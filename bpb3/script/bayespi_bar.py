import os
import argparse
import sys
import glob
from .script_high.other import common
from .script_high import parallel, bayespi_common
import numpy as np

def my_parser(parser):
#    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
#                                     description="""BayesPI-BAR delta-dbA ranking computation. This program computes """
#                                                 """rankings of given PWMs according to predicted effects of"""
#                                                 """ given sequence variants.""", add_help=False,
#                                     fromfile_prefix_chars="@")
#    try:
        required_named = parser.add_argument_group('required arguments')
        optional = parser.add_argument_group('optional arguments')

        optional.add_argument("--chemical_potentials", help="list of chemical potentials to use (e.g.,none -10 -13 -15 .. ), default is None ",
                              type=str, metavar="POTENTIAL", nargs="+")

        optional.add_argument("--reference_sequences", help="FASTA file with reference sequences, default is None",
                              type=argparse.FileType(), metavar="FASTA_FILE")

        optional.add_argument("--alternate_sequences", help="FASTA file with alternate sequences, default is None",
                              type=argparse.FileType(), metavar="FASTA_FILE")

        optional.add_argument("--pwm_folder", help="folder containing TF binding models in BayesPI .mlp format, default is None",
                              type=str, metavar="FOLDER")

        optional.add_argument("--iterations", help="iteration count for background affinity distribution calculation, default =10000 ",
                              type=int, metavar="NUMBER", default=10000)

        optional.add_argument("--p_value_cutoff", help="maximum dbA p-value to consider protein-DNA binding, default = 0.1 "
                                                       "significant",
                              type=float, metavar="NUMBER", default=0.1)

        optional.add_argument("--normalize_dba", help="divide each dbA value by its standard deviation in the"
                                                      "background sequence set, default is False",
                              action="store_true")

        optional.add_argument("--start_from_integration", help="do not compute delta-dbA values, assume they are "
                                                               "already computed. Start from integration. default is False",
                              action="store_true")

        optional.add_argument("--reuse_output", help="do not recompute dbA values if they are already present in the "
                                                     "result folder. This mode can handle partially computed results "
                                                     "from an interrupted computation. Only the missing or corrupted "
                                                     "output will be recomputed. default is False",
                              action="store_true")

        optional.add_argument("--integration", help="method to integrate delta-dbA values obtained from different"
                                                    "chemical potentials, default is pca ",
                              type=str, metavar="METHOD", default="pca", choices=["pca", "mean_ddba", "med_ddba"])

        optional.add_argument("--result_folder", help="folder to put the results (will be erased), default=result", type=str,
                              metavar="FOLDER", default="result")

        optional.add_argument("--seed", help="random seed for background affinity sampling. 0 means time-based seed, default =1 ",
                              type=int, metavar="NUMBER", default=1)

        optional.add_argument("--med_ddba_p_value", help="a threshold on \"median delta-dbA P-value\",default is None ",
                              type=float, metavar="NUMBER")

        optional.add_argument("--max_rank", help="maximum number of PWMs to write in each side of the rankings, default =10 ",
                              type=int, metavar="NUMBER", default=10)

#        optional.add_argument("-h", "--help", action="help", help="show this help message and exit")

        parallel.add_parallel_options(parser)
#
#        args = parser.parse_args()
#    except IOError as err:
#        print("Command line arguments error:", err, file=sys.stderr)
#        exit(1)
#
#    try:
#        run_bayespibar(args)
#    except RuntimeError as err:
#        print(err, file=sys.stderr)
#        exit(1)
        return parser



def get_pwms(pwm_folder):
    if not os.path.exists(pwm_folder):
        raise IOError("PWM folder, " + pwm_folder + " does not exist")

    pwm_files = glob.glob(os.path.join(pwm_folder, "*.mlp"))
    if len(pwm_files) == 0:
        raise IOError("No .mlp files found in the PWM folder " + pwm_folder)

    return pwm_files

def read_tsv_with_names(tsv_file, expected_col_names=None, expected_row_names=None):
    data, header = common.read_tsv(tsv_file, return_header=True)
    col_names = header[1:]
    if expected_col_names is not None and expected_col_names != col_names:
        raise IOError("Column names in " + tsv_file + " are not as expected")

    row_names = data[0]
    if expected_row_names is not None and expected_row_names != row_names:
        raise IOError("Row names in " + tsv_file + " are not as expected")

    data = data[1:]
    return data, col_names, row_names

def calculate_dba(fasta_file, pwm_folder, res_folder, iterations, chemical_potentials, seed, reuse_output):
    res_folder = os.path.abspath(res_folder)
    seqs = common.read_fasta(fasta_file.name)
    if len(seqs) == 0:
        raise IOError("No sequences found in " + fasta_file.name)

    seq_names = [s[0] for s in seqs]
    if len(seq_names) != len(frozenset(seq_names)):
        raise IOError("Duplicate sequence names in " + fasta_file.name + ": " +
                      ", ".join(set([x for x in seq_names if seq_names.count(x) > 1])))

    seq_names = frozenset(seq_names)

    num_sequences = len(seqs)
    pwm_files = get_pwms(pwm_folder)
    print("Calculating dbA values for", num_sequences, "sequences in", fasta_file.name, "and", len(pwm_files),
          "PWMs in", pwm_folder, "using", iterations,
          "background iterations and", len(chemical_potentials), "chemical potential(s):",
          ", ".join(chemical_potentials))

    if reuse_output and os.path.exists(res_folder):
        # print "The output folder,", res_folder, "exists. Trying to reuse the possibly partial results from it..."
        pass
    else:
        common.prepare_result_folder(res_folder)

    commands = []
    for potential in chemical_potentials:
        potential_result_folder = os.path.join(res_folder, "potential_" + potential)
        if not os.path.exists(potential_result_folder):
            os.mkdir(potential_result_folder)

        for pwm in pwm_files:
            if reuse_output:
                output_file = os.path.join(potential_result_folder, os.path.basename(pwm) +
                                           "_bindingP_0_randomBackground_" + str(iterations))

                if os.path.exists(output_file):
                    try:
                        dba_data = common.read_tsv(output_file, skip=1)
                        file_seq_names = frozenset(dba_data[0])
                        if file_seq_names == seq_names:
                            continue  # the data in the file seems ok, do not schedule the recomputation
                    except:
                        pass

            commands.append(bayespi_common.compute_affinity_command(fasta_file.name, pwm, iterations,
                                                                    potential_result_folder, potential, seed))

    return commands


def write_tsv_with_names(data_columns, id_name, col_names, row_names, file_name):
    if len(data_columns) != len(col_names):
        raise IOError("Wrong column names")

    header = [id_name] + col_names
    data_rows = [header] + list(common.transpose_list([row_names] + data_columns))
    common.write_tsv(data_rows, file_name)


def collect_dba(dba_folder):
    potentials = glob.glob(os.path.join(dba_folder, "potential_*"))

    if len(potentials) == 0:
        raise IOError("No usable dbA data found in folder " + dba_folder)

    pwms_found = None
    seqs_found = None
    for potential_folder in potentials:
        pwm_dba_files = glob.glob(os.path.join(potential_folder, "*_bindingP_0_randomBackground_*"))
        if len(pwm_dba_files) == 0:
            raise IOError("No usable dbA data found in folder " + potential_folder)

        base_pwms = frozenset(map(os.path.basename, pwm_dba_files))
        if pwms_found is None:
            pwms_found = base_pwms
        elif pwms_found != base_pwms:
            raise IOError("Different PWM results found in %s compared to %s: missing %s; extra %s" %
                          (potential_folder, potentials[0], pwms_found - base_pwms, base_pwms - pwms_found))

        p_vals = []
        dbas = []
        dba_stds = []
        pwms = []
        for pwm in pwm_dba_files:
            pwm_data = common.read_tsv(pwm, skip=1)
            if len(pwm_data) < 8:
                raise IOError("Unexpected file format for " + pwm)

            seqs = pwm_data[0]
            seqs_set = frozenset(seqs)
            if len(seqs) != len(seqs_set):
                raise IOError("Duplicate sequence names in " + pwm + ": " +
                              ", ".join(set([x for x in seqs if seqs.count(x) > 1])))

            if seqs_found is None:
                seqs_found = seqs_set
            elif seqs_found != seqs_set:
                raise IOError("Different sequence names found in " + pwm)

            p_vals.append(common.sort_by(pwm_data[4], seqs))
            dbas.append(common.sort_by(pwm_data[3], seqs))
            dba_stds.append(common.sort_by(pwm_data[7], seqs))
            pwm_name = os.path.basename(pwm)
            pwm_name = pwm_name[:pwm_name.find("_bindingP_0_randomBackground_")]
            pwms.append(pwm_name)

        p_vals = common.sort_by(p_vals, pwms)
        dbas = common.sort_by(dbas, pwms)
        dba_stds = common.sort_by(dba_stds, pwms)
        pwms = sorted(pwms)
        seqs = sorted(seqs_found)

        write_tsv_with_names(p_vals, "Sequence_name", pwms, seqs, os.path.join(potential_folder, "all_p_values.tsv"))
        write_tsv_with_names(dbas, "Sequence_name", pwms, seqs, os.path.join(potential_folder, "all_dbas.tsv"))
        write_tsv_with_names(dba_stds, "Sequence_name", pwms, seqs, os.path.join(potential_folder, "all_dba_stds.tsv"))


def get_proper_dba(dba_folder, p_value_cutoff, normalize_dba):
    p_vals, pwms, seqs = read_tsv_with_names(os.path.join(dba_folder, "all_p_values.tsv"))
    dbas, _, _ = read_tsv_with_names(os.path.join(dba_folder, "all_dbas.tsv"), pwms, seqs)
    p_vals = np.array(p_vals, float)
    dbas = np.array(dbas, float)
    dbas[dbas < 0] = 0  # negative dbA means that the TF binds to the target sequence weaker than to the background,
    # which for us is the same as no significant binding, so we replace it with 0

    dbas[p_vals >= p_value_cutoff] = 0  # the estimate of the probability that the binding to the target sequence is
    # stronger than the binding to the background must be higher than the threshold
    # otherwise, consider it as no binding

    if normalize_dba:
        dba_stds = read_tsv_with_names(os.path.join(dba_folder, "all_dba_stds.tsv"), pwms, seqs)
        dba_stds = np.array(dba_stds, float)
        std_eps = np.maximum(np.median(dba_stds, 1) * 0.1, 0.00001)
        dba_stds = np.maximum(dba_stds, std_eps[:, np.newaxis])
        dbas /= dba_stds

    return dbas, pwms, seqs


def calculate_delta_dba(ref_dba_folder, alt_dba_folder, res_folder, p_value_cutoff, normalize_dba):
    print("Calculating delta-dbA values...")
    common.prepare_result_folder(res_folder)
    potentials = list(map(os.path.basename, glob.glob(os.path.join(ref_dba_folder, "potential_*"))))
    if len(potentials) == 0:
        raise IOError("No usable dbA data found in folder " + ref_dba_folder)

    potentials_alt = map(os.path.basename, glob.glob(os.path.join(alt_dba_folder, "potential_*")))
    if frozenset(potentials) != frozenset(potentials_alt):
        raise IOError("Data for different potentials found in folders " + ref_dba_folder + " and " + alt_dba_folder)

    for potential in potentials:
        ref_potential_folder = os.path.join(ref_dba_folder, potential)
        alt_potential_folder = os.path.join(alt_dba_folder, potential)

        ref_dba, ref_pwms, ref_seqs = get_proper_dba(ref_potential_folder, p_value_cutoff, normalize_dba)
        alt_dba, alt_pwms, alt_seqs = get_proper_dba(alt_potential_folder, p_value_cutoff, normalize_dba)

        if ref_pwms != alt_pwms:
            raise IOError("PWMs differ between results in " + ref_potential_folder + " and " + alt_potential_folder)

        if ref_seqs != alt_seqs:
            raise IOError("Sequence names differ between results in " + ref_potential_folder + " and " +
                          alt_potential_folder)

        ddba = alt_dba - ref_dba  # positive change means that the affinity is increased, i.e. the alternate sequence
        # affinity is greater than the reference
        write_tsv_with_names(list(ddba), "Sequence_name", ref_pwms, ref_seqs,
                             os.path.join(res_folder, potential + "_ddba.tsv"))


def calc_mean_ddba(all_potentials_data):
    data_joined = np.array(all_potentials_data)
    return data_joined.mean(axis=0)


def calc_med_ddba(all_potentials_data):
    data_joined = np.array(all_potentials_data)
    return np.median(data_joined, axis=0)


def calc_pca_of_ddba(all_potentials_data):
    if len(all_potentials_data) == 1:
        raise IOError("Need more than one chemical potential for PCA integration")

    mean_ddba = calc_mean_ddba(all_potentials_data)
    data_joined = np.array(all_potentials_data)

    # apply a log transformation to the delta-dbA values
    log_data_joined = np.log(np.abs(data_joined) + 1) * np.sign(data_joined)

    # we process each mutation separately
    num_seqs = log_data_joined.shape[2]
    all_scores = []
    for s in range(num_seqs):
        seq_data = log_data_joined[:, :, s].transpose()
        mean_subtracted = seq_data - seq_data.mean(axis=0)[np.newaxis, :]
        U, d, Vt = np.linalg.svd(mean_subtracted, full_matrices=False)
        scores = np.dot(mean_subtracted, Vt[0])
        scores *= np.sign(np.corrcoef(mean_ddba[:, s], scores)[0, 1])
        all_scores.append(scores)

    all_scores = np.array(all_scores).transpose()
    return all_scores


def integrate_delta_dba(ddba_folder, res_folder, integration_method, potentials):
    print("Integrating delta-dbA values across chemical potentials using", integration_method, "...")

    if potentials is None:
        potentials = glob.glob(os.path.join(ddba_folder, "potential_*_ddba.tsv"))
        if len(potentials) == 0:
            raise IOError("No usable delta-dbA data found in folder " + ddba_folder)

    all_potentials_data = []
    pwms = None
    seqs = None
    pwm_output_names = None
    for potential in potentials:
        potential_ddba_file = os.path.join(ddba_folder, "potential_" + str(potential) + "_ddba.tsv")
        if not os.path.exists(potential_ddba_file):
            raise IOError("Data for chemical potential " + str(potential) + " not found in " + ddba_folder)

        ddba_data, pwms, seqs = read_tsv_with_names(potential_ddba_file, pwms, seqs)

        # remove "_disc" PWMs which are not reliable
        non_disc_pwms = [i for i, pwm in enumerate(pwms) if "_disc" not in pwm]
        ddba_data = [ddba_data[i] for i in non_disc_pwms]
        pwm_output_names = [pwms[i] for i in non_disc_pwms]

        all_potentials_data.append(np.array(ddba_data, float))

    if integration_method == "mean_ddba":
        integrated = calc_mean_ddba(all_potentials_data)
    elif integration_method == "pca":
        integrated = calc_pca_of_ddba(all_potentials_data)
    elif integration_method == "med_ddba":
        integrated = calc_med_ddba(all_potentials_data)
    else:
        raise IOError("Unknown integration method: " + integration_method)

    write_tsv_with_names(list(integrated), "Sequence_name", pwm_output_names, seqs,
                         os.path.join(res_folder, "ddba_integrated_using_" + integration_method + ".tsv"))


def descending_ranks(x):
    return np.argsort(np.argsort(-x))


def pwm_name_to_tf_tag(pwm_name):
    base = pwm_name.split("_")[0]
    if "_from_" in pwm_name:
        from_tag = pwm_name[(pwm_name.find("_from_") + 6):].split("_")[0]
    else:
        from_tag = ""

    return base + "_" + from_tag


def calculate_med_ddba_p_values(ddba_folder, pwms, seqs):
    med_ddba_file = os.path.join(ddba_folder, "ddba_integrated_using_med_ddba.tsv")
    med_ddba_data, pwms, seqs = read_tsv_with_names(med_ddba_file, pwms, seqs)
    med_ddba_data = np.array(med_ddba_data, float)

    pwm_tf_tags = map(pwm_name_to_tf_tag, pwms)
    all_nonzero_med_ddba = []
    for pwm_i in range(len(pwms)):
        for seq_i in range(len(seqs)):
            if med_ddba_data[pwm_i][seq_i] != 0:
                all_nonzero_med_ddba.append((med_ddba_data[pwm_i][seq_i], pwm_tf_tags[pwm_i], pwm_i, seq_i))

    med_ddba_p_values = np.ones_like(med_ddba_data)
    N = len(all_nonzero_med_ddba)
    if N > 1:
        all_nonzero_med_ddba.sort(key=lambda x: x[0])
        unique_tfs_up = set()
        unique_tfs_down = set()
        unique_up = np.zeros(N)
        unique_down = np.zeros(N)
        for i in range(N):
            unique_tfs_up.add(all_nonzero_med_ddba[i][1])
            unique_up[i] = len(unique_tfs_up)
            unique_tfs_down.add(all_nonzero_med_ddba[N - i - 1][1])
            unique_down[N - i - 1] = len(unique_tfs_down)

        num_positive_vals = len([x for x in all_nonzero_med_ddba if x[0] > 0])
        num_negative_vals = N - num_positive_vals
        for i in range(N):
            if all_nonzero_med_ddba[i][0] < 0:
                count = unique_up[i]
                p = float(count) / num_negative_vals
            else:
                count = unique_down[i]
                p = float(count) / num_positive_vals

            med_ddba_p_values[all_nonzero_med_ddba[i][2]][all_nonzero_med_ddba[i][3]] = p

    write_tsv_with_names(list(med_ddba_p_values), "Sequence_name", pwms, seqs,
                         os.path.join(ddba_folder, "med_ddba_p_values.tsv"))

    return med_ddba_p_values


def compute_rankings(ddba_folder, res_folder, integration_method, max_rank, med_ddba_p_value_cutoff):
    print("Computing PWM rankings for each variant...")
    common.prepare_result_folder(res_folder)
    score_ddba_file = os.path.join(ddba_folder, "ddba_integrated_using_" + integration_method + ".tsv")
    score_data, pwms, seqs = read_tsv_with_names(score_ddba_file)
    mean_ddba_file = os.path.join(ddba_folder, "ddba_integrated_using_mean_ddba.tsv")
    mean_ddba_data, _, _ = read_tsv_with_names(mean_ddba_file, pwms, seqs)

    if med_ddba_p_value_cutoff is not None:
        med_ddba_p_values = calculate_med_ddba_p_values(ddba_folder, pwms, seqs)

    pwms = np.array(pwms, str)
    mean_ddba_data = np.array(mean_ddba_data, float)
    score_data = np.array(score_data, float)

    for s in range(len(seqs)):
        seq = seqs[s]
        file_name = os.path.join(res_folder, common.normalize_string_for_file_name(seq) + "_ranking.tsv")

        def ranking_info(side):
            side_idx = np.sign(mean_ddba_data[:, s]) == np.sign(side)
            side_score = score_data[:, s][side_idx]
            side_pwms = pwms[side_idx]
            side_mean_ddba = mean_ddba_data[:, s][side_idx]
            order = np.argsort(-side_score * side)[:max_rank]
            if med_ddba_p_value_cutoff is None:
                return side_pwms[order], side_score[order], side_mean_ddba[order]

            side_med_ddba_p_value = med_ddba_p_values[:, s][side_idx]
            return side_pwms[order], side_score[order], side_mean_ddba[order], \
                   side_med_ddba_p_value[order]

        pos_info = ranking_info(1.0)
        neg_info = ranking_info(-1.0)

        def filter_by_med_ddba_p_value(arrays, mdpv):
            filtered = mdpv < med_ddba_p_value_cutoff
            return [a[filtered] for a in arrays]

        if med_ddba_p_value_cutoff is not None:
            pos_info = filter_by_med_ddba_p_value(pos_info, pos_info[3])
            neg_info = filter_by_med_ddba_p_value(neg_info, neg_info[3])

        pos_side = ["Positive"] * len(pos_info[0])
        neg_side = ["Negative"] * len(neg_info[0])
        ranking_columns = [["change"] + pos_side + neg_side,
                           ["rank"] + list(range(1, len(pos_info[0]) + 1)) + list(range(1, len(neg_info[0]) + 1)),
                           ["file_name"] + list(pos_info[0]) + list(neg_info[0]),
                           ["score_" + integration_method] + list(pos_info[1]) + list(neg_info[1]),
                           ["mean_ddba"] + list(pos_info[2]) + list(neg_info[2])]

        if med_ddba_p_value_cutoff is not None:
            ranking_columns.append(["med_ddba_p_value"] + list(pos_info[3]) + list(neg_info[3]))

        common.write_tsv(common.transpose_list(ranking_columns),
                         file_name, prefix="#\tSequence_name\t" + seq)


def run(args):
    ddba_folder = os.path.join(args.result_folder, "ddba")

    if not args.start_from_integration:
        if args.reference_sequences is None or args.alternate_sequences is None:
            raise RuntimeError("Reference and alternate sequences not specified")

        if args.chemical_potentials is None:
            raise RuntimeError("Chemical potentials not specified")

        if args.pwm_folder is None:
            raise RuntimeError("Folder with PWM files not specified")

        if args.reuse_output and os.path.exists(args.result_folder):
            print("The output folder,", args.result_folder, ", exists. Trying to reuse the (possibly partial) " +
                  "results from it...")
        else:
            common.prepare_result_folder(args.result_folder)

        ref_dba_folder = os.path.join(args.result_folder, "dba_ref")
        alt_dba_folder = os.path.join(args.result_folder, "dba_alt")

        try:
            ref_seqs = common.read_fasta(args.reference_sequences.name)
            ref_names, dups = common.uniq_and_duplicates(name for name, _ in ref_seqs)
            if len(dups) != 0:
                raise IOError("Duplicate sequence names in " + args.reference_sequences.name + ": " + ", ".join(dups))

            alt_seqs = common.read_fasta(args.alternate_sequences.name)
            alt_names, dups = common.uniq_and_duplicates(name for name, _ in alt_seqs)
            if len(dups) != 0:
                raise IOError("Duplicate sequence names in " + args.alternate_sequences.name + ": " + ", ".join(dups))

            if alt_names != ref_names:
                raise IOError("Sequence names in " + args.reference_sequences.name + " and " +
                              args.alternate_sequences.name + " are not the same")

            ref_seqs = dict(ref_seqs)
            alt_seqs = dict(alt_seqs)
            for name, ref_seq in ref_seqs.items():
                alt_seq = alt_seqs[name]
                if alt_seq.lower() == ref_seq.lower():
                    raise IOError("Same sequence in reference and alternate FASTA files: " + name)

            runner = parallel.Parallel(args, os.path.join(args.result_folder, "parallel_tmp"))
            commands = []
            commands.extend(calculate_dba(args.reference_sequences, args.pwm_folder, ref_dba_folder, args.iterations,
                                          args.chemical_potentials, args.seed, args.reuse_output))

            commands.extend(calculate_dba(args.alternate_sequences, args.pwm_folder, alt_dba_folder, args.iterations,
                                          args.chemical_potentials, args.seed, args.reuse_output))

            runner.run_commands(commands, "dba")
            print("Collecting dbA values...")
            collect_dba(ref_dba_folder)
            collect_dba(alt_dba_folder)
        except IOError as err:
            raise RuntimeError("Error computing dbA: " + err.message)

        try:
            calculate_delta_dba(ref_dba_folder, alt_dba_folder, ddba_folder, args.p_value_cutoff, args.normalize_dba)
        except IOError as err:
            raise RuntimeError("Error computing delta-dbA: " + err.args[0])

    try:
        integrate_delta_dba(ddba_folder, ddba_folder, "mean_ddba", args.chemical_potentials)
        if args.integration != "mean_ddba":
            integrate_delta_dba(ddba_folder, ddba_folder, args.integration, args.chemical_potentials)

        if args.integration != "med_ddba" and args.med_ddba_p_value is not None:
            integrate_delta_dba(ddba_folder, ddba_folder, "med_ddba", args.chemical_potentials)
    except IOError as err:
        raise RuntimeError("Error integrating delta-dbA: " + err.message)

    rankings_folder = os.path.join(args.result_folder, "rankings")
    try:
        compute_rankings(ddba_folder, rankings_folder, args.integration, args.max_rank, args.med_ddba_p_value)
    except IOError as err:
        raise RuntimeError("Error computing rankings: " + err.message)


if __name__ == "__main__":
 args=my_parser(argparse.ArgumentParser('python bayespi_bar.py')).parse_args()
 #parallel.add_parallel_options(parser)
 run(args)

