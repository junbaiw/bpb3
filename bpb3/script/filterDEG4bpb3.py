#this script is used to filter and plot differential expressed genes
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
from scipy import stats
import argparse
import os
#exec(open('filter_results4bpb3.py').read())

def my_parser(parser):
  required= parser.add_argument_group('required argument for filtering ')
  required.add_argument('--in_group1_str',help='Group1 column label', required=True, type=str, metavar="NAME")
  required.add_argument('--in_group2_str',help='Group2 column label', required=True, type=str, metavar="NAME")
  required.add_argument('--in_folder',help='Input file folder where contains bpb3 exported differential expression gene file',
                        required=True, type=str, metavar="FILE FOLDER")
  required.add_argument('--in_file', help='Input file name of bpb3 exported differential expression gene file',
                        required=True, type=str, metavar="FILE NAME")

  optional=parser.add_argument_group('optional arguments for filtering, have default values ')
  optional.add_argument('--min_median_RPKM', help='Minimum median of RPKM value in each group, default=1.0',
			 type=float, metavar="NUMBER", default=1.0)
  optional.add_argument('--rr_cutoff', help='Cutpff value for rratios to filter differential expression genes beteween two groups default=0.18,\
                        rr>0.14, 0.18, 0.22, 0.4, 0.67 represent 1.15, 1.2, 1.25, 1.5, 2.0 folder changes, respectively.',
                        type=float, metavar="NUMBER", default=0.18 )
  optional.add_argument('--sort_by_column_str', help='Column name of dataframe column that will be sorted after filtering , default=differential_expressoin_T_test_p_value',
                        type=str, metavar="NAME", default='differential_expression_T_test_p_value'  )
  optional.add_argument('--is_median', help='Use mean or median of group values to calculate folder changes between groups, default=False use mean value of the group',
                         action="store_true")
  return parser

def change_column_label(in_df,idx2column,new_label):
   tmp_cols=in_df.columns.to_list()
   tmp_cols[idx2column]=new_label
   in_df.columns=tmp_cols

def find_df_column_index(in_df,group1_str):
    col_idx=np.where(in_df.columns.str.contains(group1_str))
    return col_idx[0]

def calcult_rratio(group_mean_df,group1_str,group2_str):
   #compute relateive ratios between med-group1 and med-group2
   rratio=2*(group_mean_df[group1_str] -group_mean_df[group2_str])/(group_mean_df[group1_str] +group_mean_df[group2_str])
   group_mean_df['rratio']=list(rratio)

def calcult_group_mean_df(sorted_in_df,group1_str,group2_str,isMean):
   #calculate group mean of median of input dataframe
   #calculate group mean for each group
   group1_idx=find_df_column_index(sorted_in_df,group1_str)
   group2_idx=find_df_column_index(sorted_in_df,group2_str)

   if isMean:
     print('Use group mean for filtering ')
     group1_mean=sorted_in_df.iloc[:,group1_idx].mean(axis=1)
     group2_mean=sorted_in_df.iloc[:,group2_idx].mean(axis=1)
   else:
     print('Use group median for filtering ')
     group1_mean=sorted_in_df.iloc[:,group1_idx].median(axis=1)
     group2_mean=sorted_in_df.iloc[:,group2_idx].median(axis=1)

   group_mean_df=pd.DataFrame(columns=['gene',group1_str,group2_str])
   group_mean_df['gene']=sorted_in_df.gene
   group_mean_df[group1_str]=list(group1_mean)
   group_mean_df[group2_str]=list(group2_mean)
   return group_mean_df.copy(),group1_idx, group2_idx

def filter_df_by_rratio(merged_df, rr_cutoff,min_median_RPKM, group1_str, group2_str):
 ''' filter datafram by rratio
  filtering weak expressed genes
  rratio>0.67 -> 2 folder change
  rratio>0.4  -> 1.5 folder changes
  rratio>0.22 -> 1.25 foler changes
  rratio>0.18 -> 1.2 folder changes -->cell report
  rratio>0.14 ->  1.15 folder changes
 '''
 filtered_merged_df= merged_df[ (abs(merged_df.rratio)>=rr_cutoff) & ~((merged_df[group1_str] <min_median_RPKM ) & (merged_df[group2_str] <min_median_RPKM)) ].copy()
 filtered_merged_df=filtered_merged_df.sort_values(by='rratio')
 filtered_merged_df=filtered_merged_df.reset_index(drop=True)
 return filtered_merged_df.copy()


def zscore_of_logRPKM(filtered_merged_df, group1_idx, group2_idx,start_col_idx4sample=2 ):
  #Zscore of log transformed RPKM values for samples in datafram
  #return datafram the first column is ID and the rest of columns are data
  num_of_samples= group1_idx.shape[0]+group2_idx.shape[0]
  out4cluster_df=filtered_merged_df.iloc[:,[0]+ [i for i in range(start_col_idx4sample,1+num_of_samples+1,1)]].copy()
  norm_data=stats.zscore(np.log(out4cluster_df.iloc[:,1:].values +1.0),axis=0)
  out4cluster_df.iloc[:,1:]=norm_data
  return out4cluster_df.copy()

def run(args):
  #group1_str='HAP1_P1'
  #group2_str='HAP1_KO1'
  #in_folder='out_aug/'
  #in_file=in_folder+ group1_str+ '_vs_'+group2_str+'_differentially_expressed_genes_min1.1Fd_min1RPKM.txt'
  #min_median_RPKM=1.0
  #rr_cutoff=0.18

  group1_str=args.in_group1_str
  group2_str=args.in_group2_str
  in_folder=args.in_folder
  in_file0=args.in_file
  #below has default values
  min_median_RPKM=args.min_median_RPKM
  rr_cutoff=args.rr_cutoff
  sort_by_column_str=args.sort_by_column_str
  isMean=not args.is_median

  print('Group1 : ', group1_str, ', Group2 : ',group2_str )
  in_file=os.path.join(in_folder, in_file0)
  print('Read file : ', in_file)
  print('rr-cutoff : ', str(rr_cutoff))
  print('sort column : ', sort_by_column_str)
  print('Minimum median RPKM of group : ', str(min_median_RPKM))
  print('Is Mean of group : ', isMean)

  in_df=pd.read_csv(in_file,sep='\t')
  change_column_label(in_df,0,'gene')

  #sorted_in_df=in_df.sort_values(by=['differential_expression_T_test_p_value']).copy()
  sorted_in_df=in_df.sort_values(by=[sort_by_column_str]).copy()
  sorted_in_df=sorted_in_df.reset_index(drop=True)

  #calculate group mean for each group
  #isMean=True
  group_mean_df,group1_idx, group2_idx =calcult_group_mean_df(sorted_in_df, group1_str,group2_str,isMean)
  calcult_rratio(group_mean_df,group1_str, group2_str)
  merged_df=pd.merge(sorted_in_df,group_mean_df,on='gene',how='left')

  #filtering weak expressed genes
  filtered_merged_df=filter_df_by_rratio(merged_df, rr_cutoff,min_median_RPKM, group1_str, group2_str)

  #export filtered data
  out_file4cluster=in_file+'_rratio_filtered4cluster.csv'
  out_file4cluster=out_file4cluster.replace('.txt','')
  out_file=in_file+'_rratio_filtered.csv'
  out_file=out_file.replace('.txt','')
  print('Export at: ',out_file4cluster, out_file)

  #Zscore of log transformed RPKM values
  out4cluster_df=zscore_of_logRPKM(filtered_merged_df, group1_idx, group2_idx,start_col_idx4sample=2 )

  #export filtered data for cluster
  tmp_columns=out4cluster_df.columns.to_list()
  tmp_columns2=[ i.replace('_count','_count_zscoreLog_of') for i in tmp_columns ]
  out4cluster_df.columns=tmp_columns2
  out4cluster_df.to_csv(out_file4cluster,sep='\t',index=False)

  out_df=filtered_merged_df.copy()
  out_df.to_csv(out_file,sep='\t',index=False)
  return out_df, out4cluster_df

if __name__=='__main__':
 args=my_parser(argparse.ArgumentParser('python filterDEG4bpb3.py')).parse_args()
 run(args)



 if False:
  group1_str='HAP1_P1'
  group2_str='HAP1_KO1'
  in_folder='out_aug/' 
  in_file=in_folder+ group1_str+ '_vs_'+group2_str+'_differentially_expressed_genes_min1.1Fd_min1RPKM.txt'
  min_median_RPKM=1.0
  rr_cutoff=0.18

  in_df=pd.read_csv(in_file,sep='\t')
  change_column_label(in_df,0,'gene')

  sorted_in_df=in_df.sort_values(by=['differential_expression_T_test_p_value']).copy()
  sorted_in_df=sorted_in_df.reset_index(drop=True)

  #calculate group mean for each group
  isMean=True
  group_mean_df,group1_idx, group2_idx =calcult_group_mean_df(sorted_in_df, group1_str,group2_str,isMean)
  calcult_rratio(group_mean_df,group1_str, group2_str)
  merged_df=pd.merge(sorted_in_df,group_mean_df,on='gene',how='left')

  #filtering weak expressed genes
  filtered_merged_df=filter_df_by_rratio(merged_df, rr_cutoff,min_median_RPKM, group1_str, group2_str)

  #export filtered data
  out_file4cluster=in_file+'_rratio_filtered4cluster.csv'
  out_file4cluster=out_file4cluster.replace('.txt','')
  out_file=in_file+'_rratio_filtered.csv'
  out_file=out_file.replace('.txt','')
  print('Export at: ',out_file4cluster, out_file)

  #Zscore of log transformed RPKM values
  out4cluster_df=zscore_of_logRPKM(filtered_merged_df, group1_idx, group2_idx,start_col_idx4sample=2 )

  #export filtered data for cluster
  tmp_columns=out4cluster_df.columns.to_list()
  tmp_columns2=[ i.replace('_count','_count_zscoreLog_of') for i in tmp_columns ]
  out4cluster_df.columns=tmp_columns2
  out4cluster_df.to_csv(out_file4cluster,sep='\t',index=False)

  out_df=filtered_merged_df.copy()
  out_df.to_csv(out_file,sep='\t',index=False)





