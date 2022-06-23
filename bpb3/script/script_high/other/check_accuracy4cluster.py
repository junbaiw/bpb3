#this script is to check accuracy for cluster pwms
import glob
import os
#import check_accuracy as ck_ac
from bpb3.script.script_high.other import check_accuracy as ck_ac
import numpy as np
import pandas as pd
import argparse
#exec(open("check_accuracy4cluster.py").read())

def my_parser(parser):
   required= parser.add_argument_group("Required ")
   required.add_argument("--in_cluster_out_path", help="File path contains all clustered PWMs and uncertain PWMs, such as '../data/out/sep8_run1 ",
                          required=True, type=str, metavar="FILE PATH")
   required.add_argument("--in_bpb3_ranking_path", help="File path for bpb3 ranking results, such as ../data/out/sep8_run1/test_67snps/foreground/test_seq/rankings/",
                          required=True, type=str, metavar="FILE PATH")   
   required.add_argument("--in_fileString4clustered_pwm", help="File string or file folder/name used/contained for clustered PWMs such as quality_assessed_out_seed3",
                          required=True, type=str, metavar="FILE String")
   
   optional= parser.add_argument_group("Option parameters , have default values")
   optional.add_argument("--in_separateString4file", help="string separation for joined new file name, default= ~", 
                       type=str,metavar="string separation", default='~')
   optional.add_argument("--in_topRank_cutoff", help="Top rank cutof value for selected TFs for comparison, default=10 ",
                       type=int, metavar="NUMBER", default=10)
   optional.add_argument("--in_fileString4uncertain_pwm", help="File strin or file folder/name used/contained for untertain PWMs such as uncertain_pwms",
                        type=str, metavar="STRING NAME", default='uncertain_pwms')
   return parser

def best_rank_in_clusterResult(stf_snps,tfNameMap,record_snp2rankings,out_df,cluster_out_folder,uncertain_string,str4file_sep,quality_string):
 #for each record_snp2ranking find its mlp names for cluster
 for ri in stf_snps:
    tmp_file_name = record_snp2rankings[ri].file_name.to_list()
    record_mlps=[]

    #get tf name
    tfs, change= ck_ac.get_tfs_and_change(ri)

    #find tf name for cluster file_name
    for ii in tmp_file_name:
      if ii.startswith(quality_string):
        tmp_path=os.path.split(os.path.split(os.path.join(os.path.abspath(cluster_out_folder), ii.replace(str4file_sep,'/')))[0])[0]
        tmp_mlp=[os.path.basename(i) for i in glob.glob(os.path.join(tmp_path,'*.mlp'))]
      elif ii.startswith(uncertain_string):
        tmp_mlp=ii.split(str4file_sep)[-1]
      else:
        print('Unknonw type cluster : ' , ii )
        pass
      record_mlps.append(tmp_mlp)

    #update cluster name to tf name
    record_snp2rankings[ri]['mlps']=record_mlps
    new_row=ck_ac.find_best_rank(record_snp2rankings[ri],ri,'mlps',tfs, change,tfNameMap,out_df.columns)
    out_df=out_df.append(new_row,ignore_index=True)
 return out_df.copy()

def read_dbd_clusters(cluster_out_folder,quality_string):
 #read each DBD cluster results to a dictionary
 all_DBDs=glob.glob(os.path.join(os.path.abspath(cluster_out_folder),quality_string, '*'))
 #print(all_DBDs)
 record_DBD_clusters={}
 for di in all_DBDs:
   tmp_DBD=os.path.basename(di)
   tmp_clusters=glob.glob(os.path.join(di,'out','*'))
   record_clusters={}
   for ci in tmp_clusters:
      tmp_cluster=os.path.basename(ci)
      tmp_mlp=glob.glob(os.path.join(ci,'*.mlp'))
      record_clusters[tmp_cluster]= [os.path.basename(i) for i in tmp_mlp]
   record_DBD_clusters[tmp_DBD]=record_clusters
 return record_DBD_clusters


def run(args):
  cluster_out_folder=args.in_cluster_out_path
  bpb3_out_folder=args.in_bpb3_ranking_path
  
  quality_string=args.in_fileString4clustered_pwm
  uncertain_string=args.in_fileString4uncertain_pwm
  str4file_sep=args.in_separateString4file

  #read each DBD cluster results to a dictionary
  record_DBD_clusters=read_dbd_clusters(cluster_out_folder,quality_string)

  #read all bpb3 results 
  in_files=glob.glob(os.path.join(os.path.abspath(bpb3_out_folder),'*.tsv'))
  record_snp2rankings=ck_ac.load_results(in_files)

  #for each record_snp2ranking find its mlp names for cluster
  out_df=pd.DataFrame(columns=['mutation','bestRank','bestScore','file_name','mlps'])
  all_out_df=best_rank_in_clusterResult(list(record_snp2rankings.keys()),ck_ac.stf.tfNameMap, 
               record_snp2rankings, out_df,cluster_out_folder,uncertain_string,str4file_sep,quality_string)

  fd_out_df= best_rank_in_clusterResult(ck_ac.stf.fredrikssonMuts,ck_ac.stf.tfNameMap,
               record_snp2rankings,out_df,cluster_out_folder,uncertain_string,str4file_sep,quality_string)

  hg_out_df= best_rank_in_clusterResult(ck_ac.stf.hgmdMuts, ck_ac.stf.tfNameMap, 
               record_snp2rankings,out_df,cluster_out_folder,uncertain_string, str4file_sep,quality_string )

  ad_out_df= best_rank_in_clusterResult(ck_ac.stf.andersonEpsteinMuts,ck_ac.stf.tfNameMap, 
               record_snp2rankings,out_df,cluster_out_folder,uncertain_string,str4file_sep,quality_string)

  #calculate accuracy
  top_rank_cutoff=args.in_topRank_cutoff
  print('Considering top ', top_rank_cutoff, ' TFs , change direction: False')
  # print('AllMuts')
  # all_total_df, all_num, all_out2df=ck_ac.accuracy(all_out_df,top_rank_cutoff)

  print('FdMuts')
  fd_total_df, fd_num, fd_out2df=ck_ac.accuracy(fd_out_df,top_rank_cutoff)

  print('HgMuts')
  hg_total_df, hg_num, hg_out2df=ck_ac.accuracy(hg_out_df,top_rank_cutoff)

  print('AdMuts')
  ad_total_df, ad_num, ad_out2df=ck_ac.accuracy(ad_out_df,top_rank_cutoff)


  print('Total Muts')
  all_total_df, all_num, all_out2df=ck_ac.accuracy(all_out_df,top_rank_cutoff)


  return all_total_df, fd_total_df, hg_total_df, ad_total_df

if __name__=='__main__':
 args= my_parser(argparse.ArgumentParser('python check_accuracy4cluster.py')).parse_args()
 run(args)



