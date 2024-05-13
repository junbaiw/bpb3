import sys
import os
import argparse
from .other import common
import pandas as pd

def my_parser(parser):
#parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
#                                 description="""This program will filter those TFs whose expression is too low from """
#                                             """the results produced by affinity_change_significance_test.py. """
#                                             """Warning: it assumes a certain naming convention for PWM files. """
#                                             """If you add your own PWM files that do not conform to it it will """
#                                             """likely crash or filter them out. It is also incomplete and doesn't """
#                                             """recognize many TFs. This is the reason it is not a part of """
#                                             """the package, but rather a demo""",
#                                 add_help=False)
#try:
    required_named = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required_named.add_argument("--results_file", help="the file containing affinity_change_significance_test.py "
                                                       "output",
                                type=argparse.FileType(), required=True, metavar="FILE")

    required_named.add_argument("--expression", help="the file containing gene expression values. Two columns: "
                                                     "gene name and gene expression",
                                type=argparse.FileType(), required=True, metavar="FILE")

    optional.add_argument("--output_file", help="output file name. It has the same format as input with lowly "
                                                "expressed TFs removed", type=argparse.FileType('w'),
                          metavar="FILE", default="-")

    optional.add_argument("--min_expression", help="minimum TF gene expression value, default minimum RPKM = 0.03", type=float,
                          metavar="NUMBER", default=0.03)

#    optional.add_argument("-h", "--help", action="help", help="show this help message and exit")
#
#    args = parser.parse_args()
#except IOError as err:
#    print("Command line arguments error:", err, file=sys.stderr)
#    exit(1)
    return parser

def run(args):
 gene_rpkm = common.read_tsv(args.expression.name, skip_comments=True)
 gene_rpkm[1] = map(float, gene_rpkm[1])
 gene_rpkm[0] = [g.upper() for g in gene_rpkm[0]]
 gene_rpkm = dict(common.transpose_list(gene_rpkm))
 rpkm_threshold = args.min_expression
 
 # this map helps to infer which genes correspond to TFs from its PWM name, for the purpose of filtering out those TFs
 # whose genes are not expressed. This is specific to the naming convention for PWM files included in the demo
 # this map is incomplete! It contains those gene families/protein complexes that appear in the melanoma dataset rankings
 gene_tag_mapping = {"ETS": ["ETS1", "ETS2"], "NR1H": ["NR1H2", "NR1H3", "NR1H4"], "NFKB": ["NFKB1", "NFKB2"],
                    "NFAT": ["NFAT5"], "TATA": ["TBP"],
                    "STAT": ["STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6"],
                    "SMAD": ["SMAD" + l for l in "12345679"], "RAR": ["RARA", "RARB", "RARG"],
                    "AP1": ["FOS", "JUN"], "IRF": ["IRF" + l for l in "123456789"],
                    "AP3": ["AP3B1", "AP3B2", "AP3D1", "AP3M1", "AP3M2"], "TFAP2": ["TFAP2" + l for l in "ABCDE"],
                    "RUNX": ["RUNX" + l for l in "123"], "CDX": ["CDX" + l for l in "124"],
                    "MYF": ["MYF" + l for l in "56"], "GATA": ["GATA" + l for l in "123456"],
                    "HNF4": ["HNF4" + l for l in "AG"], "GCM": ["GCM" + l for l in "12"]}
 #add jbw
 tp_df=pd.read_csv(args.results_file.name,sep='\t')
 if tp_df.shape[0]>0:
   res_data, header = common.read_tsv(args.results_file.name, return_header=True)
   args.output_file.write("\t".join(header) + "\n")
   for i in range(len(res_data[0])):
     pwm_name = res_data[1][i]
     tag1 = pwm_name.split("_")[0]
     tag2 = pwm_name[(pwm_name.find("_from_") + 6):].split("_")[0]
     tag = tag2 if tag1 in tag2 else tag1
     if tag in gene_rpkm:
        rpkm = gene_rpkm[tag]
     elif tag1 in gene_rpkm:
        rpkm = gene_rpkm[tag1]
     elif tag2 in gene_rpkm:
        rpkm = gene_rpkm[tag2]
     elif tag in gene_tag_mapping:
        rpkm = max(gene_rpkm[t] for t in gene_tag_mapping[tag])
     elif "::" in tag:
        rpkm = max(gene_rpkm[t] for t in tag.split("::"))
     else:
        print("Unknown gene tag:", tag, ", file:", pwm_name)
        continue

     if rpkm < rpkm_threshold:
        continue

     args.output_file.write("\t".join(c[i] for c in res_data) + "\n")
 else:
   print('No restult find in : ', args.results_file.name)

if __name__=='__main__':
  args=my_parser(argparse.ArgumentParser('python filter_results_by_gene_expression.py')).parse_args()
  run(args)



