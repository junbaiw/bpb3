import argparse
import sys
import os
from .script_high import parallel
from .script_high.other import common
import numpy as np
import random
from bpb3.script import bayespi_bar
import copy
import re

def my_parser(parser):
#parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False,
#                                 fromfile_prefix_chars="@")
#try:
    optional = parser.add_argument_group('optional arguments')
    bar = parser.add_argument_group('BayesPI-BAR arguments')

    optional.add_argument("--block_size", help="the size of the block to take for background. "
                                               "It should be similar to size of the actual patient block that "
                                               "will be tested against this background, default is None",
                          type=int, metavar="NUMBER")
    #genome and regions changed file to argparse.FileType('r') by wang
    optional.add_argument("--genome", help="reference genome in FASTA format", type=argparse.FileType('r'),
                          metavar="FASTA_FILE")

    optional.add_argument("--regions", help="BED file with regions of interest from which to take the "
                                            "background samples. Regions less than --block_size  will be "
                                            "discarded. This should be similar to regions of interest that "
                                            "were used to filter the mutations", type=argparse.FileType('r'),
                          metavar="BED_FILE")

    optional.add_argument("--mutations_distribution", help="a space-separated list of numbers that represent the "
                                                           "distribution of the number of mutations to apply to "
                                                           "each background block. It should be similar to "
                                                           "distribution of mutation counts in actual patients "
                                                           "in the blocks. mussd.py outputs these distributions "
                                                           "in the last column of the block_summary.tsv file. "
                                                           "Give one number to have a constant mutation count in "
                                                           "each "
                                                           "background block, default is None", type=int, nargs="+",
                          metavar="NUMBER")

    optional.add_argument("--result_folder", help="folder to put the results (will be erased), default is result", type=str,
                          metavar="FOLDER", default="result")

    optional.add_argument("--block_resample_iterations", help="number of times to perform sampling of blocks from "
                                                              "regions. On each iteration, a number of random blocks "
                                                              "will be "
                                                              "selected from regions, default=1",
                          type=int, metavar="NUMBER", default=1)

    optional.add_argument("--block_samples_to_take", help="number of samples to take from regions on each resampling "
                                                          "iteration. For each sample, a region is chosen randomly "
                                                          "from the input regions without replacement, and then a "
                                                          "block of specified size is chosen randomly from that "
                                                          "region. If this number is not specified, it is equal to the "
                                                          "number of regions, meaning that each region will be sampled "
                                                          "once, default is None",
                          type=int, metavar="NUMBER")
    #added FileType('r') by wang
    optional.add_argument("--mutation_signature", help="file containing the mutation signature to apply when "
                                                       "generating mutations. Tab-separated two-column file with a "
                                                       "header. First column is a k-mer specification with a "
                                                       "nucleotide replacement in square brackets, i.e. A[C>A]A "
                                                       "specifies a 3-mer mutation ACA -> AAA. The second column is "
                                                       "the probability of the given replacement. Reverse-complement "
                                                       "replacement will have the same probability.",
                          type=argparse.FileType('r'), metavar="TSV_FILE")

    optional.add_argument("--background_mutations", help="file containing background mutations to use. Tab-separated "
                                                         "file without a header with at least 4 columns: chromosome, "
                                                         "position, reference nucleotide, alternate nucleotide. "
                                                         "If given, only these mutations will be selected in "
                                                         "background blocks, default is None",
                          type=argparse.FileType('r'), metavar="TSV_FILE")

#    optional.add_argument("-h", "--help", action="help", help="show this help message and exit")

    bar.add_argument("--pwm_folder", help="folder containing TF binding models in BayesPI .mlp format",
                     type=str, required=True, metavar="FOLDER")

    bar.add_argument("--chemical_potentials", help="list of chemical potentials to use",
                     type=str, metavar="POTENTIAL", required=True, nargs="+")

    bar.add_argument("--iterations", help="iteration count for background affinity distribution calculation, default = 10000",
                     type=int, metavar="NUMBER", default=10000)

    bar.add_argument("--p_value_cutoff", help="maximum dbA p-value to consider protein-DNA binding "
                                              "significant, default=0.1 ",
                     type=float, metavar="NUMBER", default=0.1)

    bar.add_argument("--normalize_dba", help="divide each dbA value by its standard deviation in the"
                                             "background sequence set, default is False",
                     action="store_true")

    bar.add_argument("--integration", help="method to integrate delta-dbA values obtained from different"
                                           "chemical potentials, default = pca",
                     type=str, metavar="METHOD", default="pca", choices=["pca", "mean_ddba", "med_ddba"])

    bar.add_argument("--seed", help="random seed for background affinity sampling. 0 means time-based seed, default is 1",
                     type=int, metavar="NUMBER", default=1)

    bar.add_argument("--start_from_integration", help="do not compute the random blocks and their delta-dbA values, "
                                                      "assume they are "
                                                      "already computed. Start from integration, default is False",
                     action="store_true")

    bar.add_argument("--reuse_output", help="do not recompute the random blocks and their dbA values if they are "
                                            "already present in the "
                                            "result folder. This mode can handle partially computed results "
                                            "from an interrupted computation. Only the missing or corrupted "
                                            "output will be recomputed, default is False",
                     action="store_true")

    bar.add_argument("--max_rank", help="maximum rank of a PWM file in the rankings to add to the background "
                                             "model , default = 10",
                          type=int, metavar="NUMBER", default=10)

    parallel.add_parallel_options(parser)

#    args = parser.parse_args()
#except IOError as err:
#    print("Command line arguments error:", err, file=sys.stderr)
#    exit(1)
    return parser

def run(args):
 #args=parser.parse_arg()

 extracted_file = os.path.join(args.result_folder, "randomly_selected_blocks.bed")
 ref_fasta = os.path.join(args.result_folder, os.path.basename(extracted_file).replace(".bed", "_ref.fasta"))
 alt_fasta = ref_fasta.replace("_ref.fasta", "_alt.fasta")

 if args.start_from_integration:
    print(
        "Skipping generation of random blocks and delta-dbA values, "
        "assuming they have already been computed and placed in",
        args.result_folder)

    if not os.path.exists(args.result_folder):
        print("The folder with precomputed blocks,", args.result_folder, "does not exist", file=sys.stderr)
        exit(1)

    if not os.path.exists(extracted_file):
        print("The file with random blocks positions,", extracted_file, "does not exist", file=sys.stderr)
        exit(1)

    if not os.path.exists(ref_fasta):
        print("The reference sequences FASTA file,", ref_fasta, "does not exist", file=sys.stderr)
        exit(1)

    if not os.path.exists(alt_fasta):
        print("The alternate sequences FASTA file,", ref_fasta, "does not exist", file=sys.stderr)
        exit(1)
 else:
    if args.genome is None:
        print("The reference genome FASTA file is not specified", file=sys.stderr)
        exit(1)

    if args.regions is None:
        print("The regions of interest BED file is not specified", file=sys.stderr)
        exit(1)

    if args.mutations_distribution is None:
        print("The distribution of mutation counts is not specified", file=sys.stderr)
        exit(1)

    if args.block_size is None:
        print("The size of background blocks is not specified", file=sys.stderr)
        exit(1)

    if args.mutation_signature and args.background_mutations:
        print("--mutation_signature and --background_mutations can not be specified at the same time", file=sys.stderr)
        exit(1)

    args.genome.close()

    if os.path.exists(args.result_folder) and args.reuse_output:
        print("Output folder", args.result_folder, ", exists. Trying to reuse data in it...")
    else:
        common.prepare_result_folder(args.result_folder)

    extract_size = args.block_size

    if not os.path.exists(extracted_file):
        print("Selecting random blocks of", extract_size, "bp as background")

        if args.background_mutations:
            print("Taking blocks from", args.regions.name, "where there is at least one background mutation")
            mutation_pos_file = os.path.join(args.result_folder, "background_mutation_pos.bed")
            os.system('bash -c "paste <(cut -f1,2 %s) <(cut -f2 %s) > %s"' %
                      (args.background_mutations.name, args.background_mutations.name, mutation_pos_file))

            genome_index_file = args.genome.name + ".fai"
            if not os.path.exists(genome_index_file):
                if not common.which("samtools"):
                    print("genome index file,", genome_index_file,
                          "does not exist. Please create it using 'samtools faidx", args.genome.name, "'",
                          file=sys.stderr)
                    exit(1)

                os.system("samtools faidx %s" % args.genome.name)

            block_file = os.path.join(args.result_folder, "possible_background_block_starts.bed")
            os.system("bedtools slop -i %s -g %s -l %d -r 0 > %s" %
                      (mutation_pos_file, genome_index_file, extract_size, block_file))

            block_file_restricted = os.path.join(args.result_folder, "possible_background_blocks_restricted.bed")
            os.system('bash -c "bedtools slop '
                      '-i <(bedtools intersect -a %s -b <(bedtools merge -i <(bedtools sort -i %s))) '
                      '-g %s -l 0 -r %d > %s"' %
                      (args.regions.name, block_file, genome_index_file, extract_size - 1, block_file_restricted))

            regions = common.read_tsv(block_file_restricted)
        else:
            regions = common.read_tsv(args.regions)

        chr_name = np.array(regions[0], str)
        start = np.array(regions[1], int)
        end = np.array(regions[2], int)
        if len(regions) > 3:
            name = regions[3]
        else:
            name = [chr_name[i] + "_" + str(start[i]) + "_" + str(end[i]) for i in range(len(chr_name))]

        if len(frozenset(name)) != len(name):
            # print len(name) - len(frozenset(name)), "duplicate gene names corrected"
            used_names = set()
            idx = 0
            for i in range(len(name)):
                if name[i] in used_names:
                    name[i] += ";dup" + str(idx)
                    idx += 1
                else:
                    used_names.add(name[i])

        name = np.array(name, str)
        upper_end = end - extract_size
        long_enough_regions = upper_end >= start
        if long_enough_regions.sum() == 0:
            print("No regions are long enough to extract a block from", file=sys.stderr)
            exit(1)

        if args.block_samples_to_take is None:
            args.block_samples_to_take = sum(long_enough_regions)

        if args.block_samples_to_take > sum(long_enough_regions):
            print("Specified --block_samples_to_take,", args.block_samples_to_take,
                  ", is greater than the number of regions long enough to contain a block,", sum(long_enough_regions),
                  file=sys.stderr)

            exit(1)

        chr_name = chr_name[long_enough_regions]
        start = start[long_enough_regions]
        end = end[long_enough_regions]
        name = name[long_enough_regions]
        upper_end = upper_end[long_enough_regions]

        start_extracted = []
        chr_extracted = []
        name_extracted = []

        if args.seed != 0:
            random.seed(args.seed)
            np.random.seed(args.seed)
        else:
            random.seed(None)
            np.random.seed(None)

        for resampling_iteration in range(args.block_resample_iterations):
            selected_regions = random.sample(range(len(start)), args.block_samples_to_take)

            for s in selected_regions:
                extracted = random.randint(start[s], upper_end[s])
                start_extracted.append(extracted)
                chr_extracted.append(chr_name[s])
                name_extracted.append(name[s] + ";+rnd" + str(extracted - start[s]))

        start_extracted = np.array(start_extracted)
        common.make_unique(name_extracted)
        common.write_tsv(
            common.transpose_list(map(list, [chr_extracted, start_extracted, start_extracted + extract_size,
                                             name_extracted])),
            extracted_file)

        print(len(name_extracted), "background blocks selected")

    if not os.path.exists(ref_fasta):
        print("Extracting reference FASTA for selected random blocks")
        common.get_fasta_by_bed(extracted_file, ref_fasta, args.genome.name)

    chr_extracted, start_extracted = common.read_tsv(extracted_file)[:2]
    chr_extracted = np.array(chr_extracted, str)
    start_extracted = np.array(start_extracted, int)


    def apply_uniform_mutations(seq, num_muts):
        for mut_pos in random.sample(range(len(seq)), num_muts):
            ref_nuc = seq[mut_pos].upper()
            if ref_nuc != "N":
                if ref_nuc not in common.nucleotide_standard_sequence:
                    print("Unknown nucleotide in reference sequence at", ref_seq[i][0], ":", ref_nuc, file=sys.stderr)
                    exit(1)

                alt_nuc = common.nucleotide_standard_sequence[(common.nucleotide_standard_sequence.index(ref_nuc) +
                                                               random.randint(1, 3)) % 4]

                seq = common.seq_replace(seq, mut_pos, alt_nuc.lower())
            else:
                seq = common.seq_replace(seq, mut_pos, common.nucleotide_standard_sequence[random.randint(0, 3)])

        return seq


    if args.mutation_signature is not None:
        muts, probs = common.read_tsv(args.mutation_signature, 1)
        src_kmers = []
        mut_kmers = []
        for mut in muts:
            src_kmers.append(re.sub(r"\[([ACGT])>[ACGT]\]", r"\1", mut.upper()))
            mut_kmers.append(re.sub(r"\[[ACGT]>([ACGT])\]", r"\1", mut.upper()))

        probs = np.array(probs, float)
        src_kmers = np.array(src_kmers, str)
        mut_kmers = np.array(mut_kmers, str)


    def find_all(a_str, sub):
        start = 0
        while True:
            start = a_str.find(sub, start)
            if start == -1:
                return

            yield start
            start += 1


    def apply_mutation_signature(seq, num_muts):
        seq = seq.upper()
        current_probs = probs
        current_src = src_kmers
        current_mut = mut_kmers
        for i in range(num_muts):
            to_keep = np.ones(len(current_probs), bool)
            for j in range(len(current_probs)):
                src = current_src[j]
                if src not in seq and common.reverse_complement(src) not in seq:
                    to_keep[j] = False

            if to_keep.sum() == 0:
                return seq

            current_probs = current_probs[to_keep]
            current_src = current_src[to_keep]
            current_mut = current_mut[to_keep]

            total_remaining_prob = current_probs.sum()
            if total_remaining_prob == 0:
                return seq

            current_probs /= total_remaining_prob

            mut_idx = np.random.choice(len(current_probs), p=current_probs)
            src = current_src[mut_idx]
            src_rc = common.reverse_complement(src)

            positions = list(find_all(seq, src)) + list(find_all(seq, src_rc))
            chosen_pos = np.random.choice(positions)
            if seq[chosen_pos:].startswith(src):
                seq = common.seq_replace(seq, chosen_pos, current_mut[mut_idx].lower())
            else:
                seq = common.seq_replace(seq, chosen_pos, common.reverse_complement(current_mut[mut_idx]).lower())

        return seq


    if args.background_mutations:
        mut_chrom, mut_pos, mut_ref, mut_alt = common.read_tsv(args.background_mutations)[:4]
        mut_chrom = np.array(mut_chrom, str)
        mut_pos = np.array(mut_pos, int)
        mut_ref = np.array(mut_ref, str)
        mut_alt = np.array(mut_alt, str)


    def apply_background_mutations(seq, chrom, pos, num_muts):
        seq = seq.upper()
        applicable_muts_idx = np.flatnonzero((mut_chrom == chrom) & (mut_pos >= pos) & (mut_pos < pos + len(seq)))
        if len(applicable_muts_idx) == 0:
            print("No applicable background mutations for %s:%d-%d" % (chrom, pos, pos + len(seq)), file=sys.stderr)
            exit(1)

        selected_muts = np.random.choice(applicable_muts_idx, min(num_muts, len(applicable_muts_idx)), replace=False)
        shifts = mut_pos[selected_muts] - pos
        m_ref = mut_ref[selected_muts]
        m_alt = mut_alt[selected_muts]
        mutated = np.zeros(len(seq), bool)
        for i in range(len(shifts)):
            if not mutated[shifts[i]] and m_ref[i].upper() != seq[shifts[i]]:
                print("Expected nucleotide %s at %s:%d, found %s (sequence pos: %d). Check that %s uses "
                      "the reference genome %s" % (m_ref[i], chrom, pos + shifts[i], seq[shifts[i]], pos,
                                                   args.background_mutations.name,
                                                   args.genome.name), file=sys.stderr)
                exit(1)

            seq = common.seq_replace(seq, shifts[i], m_alt[i].lower())
            mutated[shifts[i]] = True

        return seq


    if not os.path.exists(alt_fasta):
        ref_seq = common.read_fasta(ref_fasta)
        ref_seq = [(name, seq) for name, seq in ref_seq if common.seq_has_known(seq)]
        common.write_fasta(ref_seq, ref_fasta)
        print(len(ref_seq), "background blocks with known nucleotides")

        if args.mutation_signature:
            print("Applying supplied mutation signature to reference FASTA to make alternate FASTA")
        elif args.background_mutations:
            print("Applying supplied background mutations to reference FASTA to make alternate FASTA")
        else:
            print("Applying random uniform mutations to reference FASTA to make alternate FASTA")

        for i in range(len(ref_seq)):
            num_muts = random.sample(args.mutations_distribution, 1)[0]
            seq = ref_seq[i][1]

            if args.mutation_signature:
                mutated_seq = apply_mutation_signature(seq, num_muts)
            #elif args.background_mutations:
                mutated_seq = apply_background_mutations(seq, chr_extracted[i], start_extracted[i] + 1, num_muts)
            else:
                mutated_seq = apply_uniform_mutations(seq, num_muts)

            ref_seq[i] = ref_seq[i][0], mutated_seq

        common.write_fasta(ref_seq, alt_fasta)

 print("Running BayesPI-BAR to compute background affinity changes. This will take a while...")
 #changed by wang
 #bayespi_bar_args = copy.deepcopy(args)
 #bayespi_bar_args=parser.parse_args()
 bayespi_bar_args=argparse.Namespace(**vars(args))
 for attr, value in vars(args).items():
    try:
        setattr(bayespi_bar_args, attr, copy.deepcopy(value))
    except Exception:
        #pass
        # Alternatively, copy reference when copy impossible:
        print(attr,value)
        setattr(bayespi_bar_args, attr, value)
 #end wang
 setattr(bayespi_bar_args, "result_folder", os.path.join(args.result_folder, "bayespi_bar_result"))
 setattr(bayespi_bar_args, "reference_sequences", open(ref_fasta))
 setattr(bayespi_bar_args, "alternate_sequences", open(alt_fasta))
 setattr(bayespi_bar_args, "med_ddba_p_value", None)
 #setattr(bayespi_bar_args, "max_rank", 10)
 bayespi_bar.run(bayespi_bar_args)

if __name__=='__main__':
 args=my_parser(argparse.ArgumentParser('python background_affinity_changes.py')).parse_args()
 run(args)






