import glob
import os
from .script_high.other import common
import collections
import argparse
import pathlib 
import pandas as pd

#melanoma_folder=os.path.abspath("../../data/skin_cancer/out")
#foreground_folder = os.path.join(melanoma_folder, "foreground")
#background_folder = os.path.join(melanoma_folder, "background")
#mussd_result_folder = os.path.join(melanoma_folder, "mussd_blocks")

def my_parser(parser):
  required= parser.add_argument_group('Required')
  required.add_argument('--foreground_folder', help='Folder path for result of foreground ',
                         required=True, type=pathlib.Path, metavar="FOLDER")
  required.add_argument('--background_folder', help="Folder path for result of background ",
                         required=True, type=pathlib.Path, metavar="FOLDER")
  required.add_argument('--mussd_result_folder', help="Folder path for result of mussd ",
                         required=True, type=pathlib.Path, metavar="FOLDER")
  return parser

def pwm_gene_tag(pwm_name):
    tag1 = pwm_name.split("_")[0]
    tag2 = pwm_name[(pwm_name.find("_from_") + 6):].split("_")[0]
    tag = tag2 if "," not in tag2 and len(tag2) < 10 else tag1
    return tag.upper()

def run(args):
 """Here is default seting of files under foregroud folder which will be used to generate figures """
 foreground_folder=os.path.abspath(args.foreground_folder)
 background_folder=os.path.abspath(args.background_folder)
 mussd_result_folder=os.path.abspath(args.mussd_result_folder)

 blocks_summary = common.read_tsv(os.path.join(mussd_result_folder, "blocks_summary.tsv"), skip=1)
 blocks_summary = dict(zip(blocks_summary[0], common.transpose_list(blocks_summary[1:])))

 for res_folder in glob.glob(os.path.join(foreground_folder, "*")):
   print(res_folder)
   #add jbw
   res_file = os.path.join(res_folder, "result_filtered.tsv")
   #tmp_pd=pd.read_csv(res_file,sep='\t')
   #if os.stat(res_file).st_size>0 and tmp_pd.shape[0]>0:
   if common.is_non_zero_file(res_file) and os.stat(res_file).st_size>0:
     tmp_pd=pd.read_csv(res_file,sep='\t')
     if tmp_pd.shape[0]>0:
       columns, header = common.read_tsv(res_file, return_header=True)

       pwm_names = columns[1]
       gene_names = list(map(pwm_gene_tag, pwm_names))
       gene_name_counts = collections.Counter(gene_names)
       for name, count in gene_name_counts.items():
         if count > 1:
            idx = 1
            for i in range(len(gene_names)):
                if gene_names[i] == name:
                    gene_names[i] += " (" + str(idx) + ")"
                    idx += 1

       columns.append(gene_names)
       header.append("gene_name")
       renamed_file = os.path.join(res_folder, "result_gene_names.tsv")
       common.write_tsv(common.transpose_list(columns), renamed_file, prefix="\t".join(header))

       block_id = os.path.basename(res_folder)
       block_data = blocks_summary[block_id]
       title = "chr" + block_data[0] + ":" + block_data[1] + "-" + block_data[2]
       if len(block_data) > 6:
         region_names = ", ".join(r.split(";")[0] for r in block_data[6].split(","))
         title += " " + region_names

       result_code=os.system("bpb3 plot_result " +
              " --background_folder " + background_folder +
              " --foreground_folder " + res_folder +
              " --significant_pwms " + renamed_file +
              " --title " + common.quote(title) +
              " --output_file " + os.path.join(res_folder, "result_plot.png") +
              " --use_tags")
       if result_code !=0:
         print("Plot file error in ", res_folder)
         exit(1)
     else:
         print('No result find in : result_filtered.rsv' )
   else:
     print('No result find in : result_filtered.rsv' )

if __name__=='__main__':
  args=my_parser(argparse.ArgumentParser('python make_plots.py')).parse_args()
  run(args)



