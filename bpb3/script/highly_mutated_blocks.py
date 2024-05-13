import argparse
import sys
from .script_high.other import common
import os
import scipy.stats
import numpy as np
import pandas as pd

def my_parser(parser):
#parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
#                                 description="""This program finds which blocks produced by mussd.py have """
#                                             """significantly more mutations than would be expected if all mutations """
#                                             """were uniformly distributed across regions of interest""",
#                                 add_help=False)
#try:
    required_named = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required_named.add_argument("--blocks_folder", help="the folder containing mussd.py output",
                                type=str, required=True, metavar="FOLDER")

    optional.add_argument("--output_file", help="output file name. Each line will contain an ID of a significantly "
                                                "mutated block, default is -", type=argparse.FileType('w'),
                          metavar="FILE", default="-")

    optional.add_argument("--p_value", help="P-value cutoff for the binomial test, default=0.05", type=float,
                          metavar="NUMBER", default=0.05)
    #added jbw
    optional.add_argument("--correctPval", help="Whether adjust P-values by bonferroni correction, default is True, if use False then this option will be turned off for correction",
                         action="store_true")
    optional.add_argument("--useRegionCounts", help="Whether use block counts or region counts to calculate expected probability that mutation falls in a certain region. Default=block counts", action="store_true")
    optional.add_argument("--useMutationsInRegions",help="Whether use mutations in regions or mutations in blocks to do bionormal test. Default=mutations in bĺocks", action="store_true")
    #optional.add_argument("-h", "--help", action="help", help="show this help message and exit")


#    args = parser.parse_args()
#except IOError as err:
#    print("Command line arguments error:", err, file=sys.stderr)
#    exit(1)
    return parser

def run(args):
 #added jbw
 is_correctPval=args.correctPval

 if not os.path.exists(args.blocks_folder):
    print("Folder does not exist:", args.blocks_folder, file=sys.stderr)
    exit(1)

 muts_file = os.path.join(args.blocks_folder, "mutations_summary.tsv")
 blocks_file = os.path.join(args.blocks_folder, "blocks_summary.tsv")
 if not os.path.exists(muts_file) or not os.path.exists(blocks_file):
    print("Folder", args.blocks_folder, "does not contain proper mussd.py output", file=sys.stderr)
    exit(1)

 info = dict(common.transpose_list(common.read_tsv(muts_file)))
 #added by wang 2021 feb
 if ('regions_count' not in info.keys()) | (not args.useRegionCounts):
   print('no regions_count but use blocks_count instead!')
   regions_count=int(info['blocks_count'])
 else:
   regions_count = int(info["regions_count"])

  #added jbw 2024
 if ('mutations_in_regions' not in info.keys()) | (not args.useMutationsInRegions):
   print('no mutations_in_regions but use mutations_in_blocks instead!')
   mutations_in_regions= int(info['mutations_in_blocks'])
 else:
   mutations_in_regions = int(info["mutations_in_regions"])

 #add jbw
 if os.stat(blocks_file).st_size>0 and pd.read_csv(blocks_file,sep='\t').shape[0]>0 :
   print('Read ,', blocks_file)
   blocks_data = common.read_tsv(blocks_file, skip=1)
 else:
   print("No blocks find in , ", blocks_file, ' program exit' )
   exit(1)
 num_muts = map(int, blocks_data[4])
 list2num_muts=list(num_muts)
 p_vals_corrected=[]
 expected_probability_that_mutation_falls_in_a_certain_region = 1.0 / regions_count
 #added jbw
 if is_correctPval:
    print('Adjust P-values by Bonferroni correction in high mutation blocks test')
    #print( mutations_in_regions, expected_probability_that_mutation_falls_in_a_certain_region, regions_count)
    p_vals_corrected = [scipy.stats.binom_test(ni, mutations_in_regions,  expected_probability_that_mutation_falls_in_a_certain_region)*regions_count for ni in list2num_muts]
    np_val=np.array(p_vals_corrected)
    np_val[np_val>1]=1
    p_vals_corrected=list(np_val)
 else:
    print('No P-value correction is applied in high mutation blocks test')
    p_vals_corrected = [scipy.stats.binom_test(ni, mutations_in_regions, expected_probability_that_mutation_falls_in_a_certain_region)  for ni in list2num_muts]
    

 for i in range(len(p_vals_corrected)):
    if p_vals_corrected[i] < args.p_value:
        args.output_file.write(blocks_data[0][i] + "\n")

if __name__=='__main__':
  args=my_parser(argparse.ArgumentParser('python highly_mutated_blocks.py')).parse_args()
  run(args)


