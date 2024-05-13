import os
import sys
import argparse
from .script_high.other import common
import collections
import glob
import shutil
import numpy as np
from .script_high import parallel

def my_parser(parser):
#parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
#                                 description="""This program selects parameters (block size, mutation distribution, """
#                                             """set of PWM files) """
#                                             """for background computation """
#                                             """based on selected mutation blocks obtained by mussd.py. """
#                                             """It will create a shell script with the command to create the """
#                                             """background model with selected parameters.""",
#                                 add_help=False,
#                                 fromfile_prefix_chars="@")
#try:
    #added FileType('r') by wang
    required_named = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    bar = parser.add_argument_group('BayesPI-BAR arguments')

    required_named.add_argument("--blocks_folder", help="the folder containing mussd.py output",
                                type=str, required=True, metavar="FOLDER")

    required_named.add_argument("--block_ids", help="the space-separated list of block IDs against which the "
                                                    "background model is being made",
                                type=str, required=True, metavar="ID", nargs="+")

    required_named.add_argument("--foreground_folder", help="the folder containing bayespi_bar.py output for each "
                                                            "selected block",
                                type=str, required=True, metavar="FOLDER")

    required_named.add_argument("--pwm_folder", help="folder containing TF binding models in BayesPI .mlp format",
                                type=str, required=True, metavar="FOLDER")

    required_named.add_argument("--background_pwm_folder", help="folder where to put PWM files selected for background",
                                type=str, required=True, metavar="FOLDER")

    required_named.add_argument("--background_output_folder", help="folder where to put the background model",
                                type=str, required=True, metavar="FOLDER")

    required_named.add_argument("--genome", help="reference genome in FASTA format", type=argparse.FileType('r'),
                                required=True,
                                metavar="FASTA_FILE")

    required_named.add_argument("--regions", help="BED file with regions of interest from which to take the "
                                                  "background samples. This should be similar to regions of interest "
                                                  "that "
                                                  "were used to filter the mutations", type=argparse.FileType('r'),
                                required=True,
                                metavar="BED_FILE")

    optional.add_argument("--output_file", help="output file name. It will contain the command to create the "
                                                "background model, default is -", type=argparse.FileType('w'),
                          metavar="FILE", default="-")

    optional.add_argument("--max_rank", help="maximum rank of a PWM file in the rankings to add to the background "
                                             "model , default = 30",
                          type=int, metavar="NUMBER", default=30)

    optional.add_argument("--block_resample_iterations", help="number of times to perform sampling of blocks from "
                                                              "regions. On each iteration, a number of random blocks "
                                                              "will be "
                                                              "selected from regions, default =1 ",
                          type=int, metavar="NUMBER", default=1)

    optional.add_argument("--block_samples_to_take", help="number of samples to take from regions on each resampling "
                                                          "iteration. For each sample, a region is chosen randomly "
                                                          "from the input regions without replacement, and then a "
                                                          "block of specified size is chosen randomly from that "
                                                          "region. If this number is not specified, it is equal to the "
                                                          "number of regions, meaning that each region will be sampled "
                                                          "once, default is None",
                          type=int, metavar="NUMBER")

    optional.add_argument("--mutation_signature", help="file containing the mutation signature to apply when "
                                                       "generating mutations. Tab-separated two-column file with a "
                                                       "header. First column is a k-mer specification with a "
                                                       "nucleotide replacement in square brackets, i.e. A[C>A]A "
                                                       "specifies a 3-mer mutation ACA -> AAA. The second column is "
                                                       "the probability of the given replacement. Reverse-complement "
                                                       "replacement will have the same probability. Default is None ",
                          type=argparse.FileType('r'), metavar="TSV_FILE")

    optional.add_argument("--background_mutations", help="file containing background mutations to use. Tab-separated "
                                                         "file without a header with at least 4 columns: chromosome, "
                                                         "position, reference nucleotide, alternate nucleotide. "
                                                         "If given, only these mutations will be selected in "
                                                         "background blocks, default is None",
                          type=argparse.FileType('r'), metavar="TSV_FILE")

    optional.add_argument("--isCluster4pwm", help="whether input PWMs are clustered PWMs or not, default is False , original PWMs without cluster",
                           action="store_true")

    #optional.add_argument("-h", "--help", action="help", help="show this help message and exit")

    bar.add_argument("--chemical_potentials", help="list of chemical potentials to use",
                     type=str, metavar="POTENTIAL", required=True, nargs="+")

    bar.add_argument("--iterations", help="iteration count for background affinity distribution calculation, default=10000",
                     type=int, metavar="NUMBER",default=10000)

    bar.add_argument("--p_value_cutoff", help="maximum dbA p-value to consider protein-DNA binding significanti, default=0.1",
                     type=float, metavar="NUMBER", default=0.1)

    bar.add_argument("--reuse_output", help="do not recompute the random blocks and their dbA values if they are "
                                            "already present in the "
                                            "result folder. This mode can handle partially computed results "
                                            "from an interrupted computation. Only the missing or corrupted "
                                            "output will be recomputed, default is False",
                     action="store_true")
    
    bar.add_argument("--start_from_integration", help="do not compute the random blocks and their delta-dbA values, "
                                                      "assume they are "
                                                      "already computed. Start from integration, default is False",
                     action="store_true")

    parallel.add_parallel_options(parser)

#    args = parser.parse_args()
#except IOError as err:
#    print("Command line arguments error:", err, file=sys.stderr)
#    exit(1)
    return parser 

def run(args):

 if not os.path.exists(args.blocks_folder):
    print("Folder does not exist:", args.blocks_folder, file=sys.stderr)
    exit(1)

 blocks_file = os.path.join(args.blocks_folder, "blocks_summary.tsv")
 if not os.path.exists(blocks_file):
    print("Folder", args.blocks_folder, "does not contain proper mussd.py output", file=sys.stderr)
    exit(1)

 #added wang
 args.genome.close()
 args.regions.close()

 significant_block_tags = frozenset(args.block_ids)
 significant_block_sizes = []
 significant_block_mutations_distribution = []
 blocks_data = common.read_tsv(blocks_file, skip=1)
 for i in range(len(blocks_data[0])):
    if blocks_data[0][i] in significant_block_tags:
        significant_block_sizes.append(int(blocks_data[3][i]) - int(blocks_data[2][i]))
        significant_block_mutations_distribution.extend(map(int, blocks_data[6][i].split(",")))

 if len(significant_block_sizes) == 0:
    print("Specified blocks not found in", blocks_file, file=sys.stderr)
    exit(1)

 med_block_size = int(np.median(significant_block_sizes))
 print("Background block size will be equal to median significant block size:", med_block_size)
 print("Mutation count distribution:")
 counts = collections.Counter(significant_block_mutations_distribution)
 print("\t", "number of mutations", "\t", "probability")
 total_count = sum(counts.values())
 for n in sorted(counts.keys()):
    print("\t", n, "\t", float(counts[n]) / total_count)

 pwms_for_background = set()

 for block_tag in significant_block_tags:
    ranking_folder = os.path.join(args.foreground_folder, block_tag, "rankings")
    for ranking_file in glob.glob(os.path.join(ranking_folder, "*.tsv")):
        data = common.read_tsv(ranking_file, skip=2)
        rank = np.array(data[1], int)
        pwm = np.array(data[2], str)
        pwms_for_background.update(pwm[rank <= args.max_rank])

 print(len(pwms_for_background), "PWMs are chosen for background model based on foreground rankings (max rank is",
      args.max_rank, ")")
 common.clear_folder(args.background_pwm_folder)
 print(args.background_pwm_folder,args.pwm_folder)
 for pwm in pwms_for_background:
    if args.isCluster4pwm:
      #change file name before copy (remove cluster and DBD name from the file name)
      pwm1=pwm.split('_')
      if pwm1[0].isnumeric():
        pwm1='_'.join(pwm1[1:])
      else:
        pwm1=pwm 
      pwm2='-'.join(pwm1.split('-')[0:-1])+'.mlp'
      shutil.copy(os.path.join(args.pwm_folder, pwm2), args.background_pwm_folder)
    else:
      shutil.copy(os.path.join(args.pwm_folder, pwm), args.background_pwm_folder)

 sbatch_header = parallel.make_sbatch_header("bg_master", args.slurm_account, 200, 1500, 1)
 #removed name in genome, regions, mutation_signature, background_mutations by wang
 #if args.start_from_integration:
 #   str2start_from_integration=" --start_from_integration " 
 #else:
 #   str2start_from_integration=""

 background_command = "bpb3 background_affinity_changes " + " \\\n\t --block_size " + str(
    med_block_size) + " \\\n\t" \
                      " --result_folder " + args.background_output_folder + " \\\n\t" + \
                     " --genome " + args.genome.name + " \\\n\t" + \
                     " --regions " + args.regions.name + " \\\n\t" + \
                     " --mutations_distribution " + " ".join(map(str, significant_block_mutations_distribution)) + \
                     " \\\n\t" + \
                     " --max_rank " + str(args.max_rank) + " \\\n\t" + \
                     " --start_from_integration" * args.start_from_integration + " \\\n\t" + \
                     " --pwm_folder " + args.background_pwm_folder + " \\\n\t" + \
                     " --chemical_potentials " + " ".join(args.chemical_potentials) + " \\\n\t" + \
                     ((" --iterations " + str(args.iterations)) + " \\\n\t" if args.iterations is not None else "") + \
                     ((" --p_value_cutoff " + str(args.p_value_cutoff)) + " \\\n\t"
                      if args.p_value_cutoff is not None else "") + " \\\n\t" + \
                     " --integration mean_ddba" + " \\\n\t" + \
                     " --block_resample_iterations " + str(args.block_resample_iterations) + " \\\n\t" + \
                     ((" --block_samples_to_take " + str(args.block_samples_to_take))
                      if args.block_samples_to_take is not None else "") + " \\\n\t" + \
                     " --use_cores " + str(args.use_cores) + " \\\n\t" + \
                     (" --use_slurm \\\n\t" if args.use_slurm else "") + \
                     ((" --slurm_account " + args.slurm_account) if args.slurm_account is not None else "") + \
                     " --max_nodes " + str(args.max_nodes) + " \\\n\t" + \
                     (" --reuse_output \\\n\t" if args.reuse_output else "") + \
                     ((" --mutation_signature %s" % args.mutation_signature.name)
                      if args.mutation_signature is not None else "") + \
                     ((" --background_mutations %s" % args.background_mutations.name)
                      if args.background_mutations is not None else "")

 args.output_file.write(sbatch_header + "\n" + background_command + "\n")


if __name__ =='__main__':
 args=my_parser(argparse.ArgumentParser('python choose_background_parameters.py')).parse_args()
 run(args)
