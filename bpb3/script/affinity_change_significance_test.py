import os
from .script_high.other import common
import glob
import numpy as np
import argparse
import sys
from .script_high.other.test import mannwhitneyu

def my_parser(parser):
#    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)
#    try:
        required_named = parser.add_argument_group('required arguments')
        optional = parser.add_argument_group('optional arguments')

        required_named.add_argument("--background_folder", help="folder containing the computed background model "
                                                                "produced by background_affinity_changes.py",
                                    type=str, required=True, metavar="FOLDER")

        required_named.add_argument("--foreground_folder", help="folder containing the computed delta-dba values for "
                                                                "patient data produced by bayespi_bar.py",
                                    type=str, required=True, metavar="FOLDER")

        optional.add_argument("--output_file", help="output file name, default is -", type=argparse.FileType('w', encoding="utf-8"),
                              metavar="FILE", default="-")

        optional.add_argument("--p_value", help="P-value cutoff for the Wilcoxon test, default=0.05", type=float,
                              metavar="NUMBER", default=0.05)

        optional.add_argument("--max_rank", help="maximum rank of PWM to consider in the patient results, default=15", type=int,
                              metavar="NUMBER", default=15)

        optional.add_argument("--exact_test",
                              help="use exact Wilcoxon rank-sum test from R. R needs to be installed, default is False",
                              action="store_true")
        optional.add_argument("--pval_correction", help="Whether adjust P-values by bonferroni correction, default is False",
                              action="store_true")

#        optional.add_argument("-h", "--help", action="help", help="show this help message and exit")
#
#        args = parser.parse_args()
#    except IOError as err:
#        print("Command line arguments error:", err, file=sys.stderr)
#        exit(1)
        return parser


def load_mean_ddba(ddba_file_name):
    mean_ddba_data, header = common.read_tsv(ddba_file_name, return_header=True)
    mean_ddba_data = np.array(mean_ddba_data[1:], float)
    header = header[1:]

    mean_ddba_distribution = {}
    for file_i in range(mean_ddba_data.shape[0]):
        file_name = header[file_i]
        if "_disc" in file_name:
            continue

        mean_ddba_distribution[file_name] = mean_ddba_data[file_i, :]

    print(ddba_file_name, ":", mean_ddba_data.shape[1], "mutations, ", mean_ddba_data.shape[0], "PWMs")
    return mean_ddba_distribution


def pwm_tf_tag(pwm_file_name):
    fields = pwm_file_name.split("_")
    n1 = fields[0].lower()
    if "from" in fields:
        n2 = fields[fields.index("from") + 1].lower()
    else:
        n2 = n1

    return n1 + ":" + n2


def get_posneg_files_from_rankings(ranking_folder, max_rank):
    positive_files = set()
    negative_files = set()

    for mut_ranking_file in glob.glob(os.path.join(ranking_folder, "*.tsv")):
        ranks = common.transpose_list(common.read_tsv(mut_ranking_file, skip_comments=True, skip=1))
        # print os.path.basename(ranks_file)
        for pwm_info in ranks:
            change, rank, file_name = pwm_info[:3]
            if max_rank is not None and int(rank) > max_rank:
                continue

            (positive_files if change == "Positive" else negative_files).add(file_name)

    return positive_files, negative_files


def compute_genome_wide_significance(foreground_result_folder, background_result_folder, out_file, p_cutoff, max_rank,
                                     use_exact_test, pval_correction):
    if pval_correction:
      out_file.write("direction\tpwm_name\tbonferroni_p\taffected_cases\ttotal_cases\n")
    else:
      out_file.write("direction\tpwm_name\tuncorrected_p\taffected_cases\ttotal_cases\n")      

    background_mean_ddba_file = os.path.join(background_result_folder, "ddba", "ddba_integrated_using_mean_ddba.tsv")
    mean_ddba_distribution_background = load_mean_ddba(background_mean_ddba_file)

    foreground_mean_ddba_file = os.path.join(foreground_result_folder, "ddba", "ddba_integrated_using_mean_ddba.tsv")
    mean_ddba_distribution_real = load_mean_ddba(foreground_mean_ddba_file)
    num_pwms = len(mean_ddba_distribution_real)
    print(num_pwms, "total PWMs")

    foreground_ranking_folder = os.path.join(foreground_result_folder, "rankings")
    positive_files, negative_files = get_posneg_files_from_rankings(foreground_ranking_folder, max_rank)

    significant_positive_files = []
    not_found = set()
    for file_name in positive_files:
        if file_name in mean_ddba_distribution_real:
            if file_name not in mean_ddba_distribution_background:
                not_found.add(file_name)
                # print "Background data for", file_name, "was not computed. Skipping significance test."
                continue

            if use_exact_test:
                genome_p = common.exact_wilcox_test(mean_ddba_distribution_background[file_name],
                                                    mean_ddba_distribution_real[file_name],
                                                    side="less")
                if pval_correction:
                  genome_p_corrected = genome_p * num_pwms
                else:
                  genome_p_corrected = genome_p
            else:
                genome_p = mannwhitneyu(mean_ddba_distribution_background[file_name],
                                        mean_ddba_distribution_real[file_name],
                                        alternative="less")
                if pval_correction:
                  genome_p_corrected = genome_p.pvalue * num_pwms
                else:
                  genome_p_corrected = genome_p.pvalue

            if p_cutoff < 0 or genome_p_corrected < p_cutoff:
                significant_positive_files.append((file_name, genome_p_corrected,
                                                   sum(mean_ddba_distribution_real[file_name] > 0)))
                # b = mean_ddba_distribution_background[file_name]
                # print file_name, ":", ", ".join(map(str, b[b!=0]))

    significant_positive_files.sort(key=lambda x: x[1])
    if pval_correction :
       print("Positive PWMs at", p_cutoff, "bonferroni significance:")
    else:
       print("Positive PWMs at", p_cutoff, "uncorrected pval significance:")

    for file_name, genome_p_corrected, num_cases in significant_positive_files:
        if pval_correction:
          print("\t", file_name, "\tbonferroni p:", genome_p_corrected, "\t#cases:", num_cases, "/",
               len(mean_ddba_distribution_real[file_name]))
        else:
          print("\t", file_name, "\tuncorrected p:", genome_p_corrected, "\t#cases:", num_cases, "/",
               len(mean_ddba_distribution_real[file_name]))

        out_file.write(
            ("Positive\t" + file_name + "\t" +
             str(genome_p_corrected) +
             "\t" + str(num_cases) + "\t" + str(len(mean_ddba_distribution_real[file_name])) +
             "\n"))

    significant_negative_files = []
    for file_name in negative_files:
        if file_name in mean_ddba_distribution_real:
            if file_name not in mean_ddba_distribution_background:
                not_found.add(file_name)
                # print "Background data for", file_name, "was not computed. Skipping significance test."
                continue

            if use_exact_test:
                genome_p = common.exact_wilcox_test(mean_ddba_distribution_background[file_name],
                                                    mean_ddba_distribution_real[file_name],
                                                    side="greater")
                if pval_correction:
                   genome_p_corrected = genome_p * num_pwms
                else:
                   genome_p_corrected = genome_p 
                  
            else:
                genome_p = mannwhitneyu(mean_ddba_distribution_background[file_name],
                                        mean_ddba_distribution_real[file_name],
                                        alternative="greater")
              
                if pval_correction:
                   genome_p_corrected = genome_p.pvalue * num_pwms
                else:
                   genome_p_corrected = genome_p.pvalue

            if p_cutoff < 0 or genome_p_corrected < p_cutoff:
                significant_negative_files.append((file_name, genome_p_corrected,
                                                   sum(mean_ddba_distribution_real[file_name] < 0)))

        # print file_name, genome_p_corrected, ":", ", ".join(map(str, mean_ddba_distribution_real[file_name]))

    significant_negative_files.sort(key=lambda x: x[1])
    if pval_correction:
       print("Negative PWMs at", p_cutoff, "bonferroni significance:")
    else:
       print("Negative PWMs at", p_cutoff, "uncorrected pval significance:")

    for file_name, genome_p_corrected, num_cases in significant_negative_files:
        if pval_correction:
          print("\t", file_name, "\tbonferroni p:", genome_p_corrected, "\t#cases:", num_cases, "/",
              len(mean_ddba_distribution_real[file_name]))
        else:
          print("\t", file_name, "\tuncorrected p:", genome_p_corrected, "\t#cases:", num_cases, "/",
              len(mean_ddba_distribution_real[file_name]))

        
        out_file.write(
            ("Negative\t" + file_name + "\t" +
             str(genome_p_corrected) +
             "\t" + str(num_cases) + "\t" + str(len(mean_ddba_distribution_real[file_name])) +
             "\n")) #.encode("utf-8"))

    if len(not_found) > 0:
        print("For the following PWMs, the background data was not found, so they were not tested:",
              ", ".join(sorted(not_found)))

def run(args):
    compute_genome_wide_significance(args.foreground_folder, args.background_folder, args.output_file, args.p_value,
                                     args.max_rank, args.exact_test, args.pval_correction)


if __name__ == "__main__":
  args=my_parser(argparse.ArgumentParser('python affinity_change_significance_test.py ')).parse_args()
  run(args)

