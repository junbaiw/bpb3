#this script is used to compute mutation scores for SNPs that evaluated by bpb3 with top ranked TFs
import pandas as pd
import glob
import os
import argparse

def my_parser(parser):
  required_named= parser.add_argument_group('required arguments')
  optional= parser.add_argument_group('optional arguments')

  required_named.add_argument('-in_folder','--in_folder_path',help='input TF ranked mutation blocks/SNPs path that exported by bpb3 package', 
                              type=str, required=True, metavar='IN_FOLDER_PATH')
  optional.add_argument('-in_file_str','--in_file_postprefix_string', help='input file postprefix string , default=*.tsv',
                              type=str, default='*.tsv')
  optional.add_argument('-rank','--top_rank_TFs',  
      help='the top N ranked TFs (both positive and negative changes) that will be used to compute mutation scores (mean of absolute scores), default N=20',
                         type=int, default=20)
  optional.add_argument('-score','--select_bpb3_score', help='select bpb3 exported scores for computing mutation scores, default=score_pca ',
                         type=str, default='score_pca')

  return parser

def run(args):
 #in_folder='rankings/bpb3'
 #in_file_postprefix='*.tsv'
 #top_rank=20
 
 in_folder=args.in_folder_path
 in_file_postprefix=args.in_file_postprefix_string
 top_rank=args.top_rank_TFs
 select_score=args.select_bpb3_score

 in_files=glob.glob(os.path.join(in_folder, in_file_postprefix))
 record_ki=[]
 record_data=[]
 #loop in each SNP export file
 loop=0
 print('Use top ', str(top_rank), ' TFs ', select_score, ' to compute mutation scores!')
 for fi in in_files:
   #print(fi)
   tmp_df=pd.read_csv(fi,sep='\t',skiprows=1)
   tmp_data=tmp_df.loc[tmp_df['rank']<=top_rank,select_score].abs().mean()
   #print(tmp_data)
   #record snp or mutation block file name
   record_ki.append(os.path.basename(fi))
   record_data.append(tmp_data)
   loop +=1
  
 print('Total ', str(loop), ' mutation blocks/SNPs loaded.')
 #export mutation scores for SNPs or mutation blocks
 out_df=pd.DataFrame(data=[record_ki,record_data]).T
 out_df.columns=['id','score']
 out_file=os.path.join(in_folder,'allMutation.scores')
 print(out_file)
 out_df.to_csv(out_file,sep='\t',index=False)


if __name__=='__main__':
  args=my_parser(argparse.ArgumentParser('python compute_mutation_score.py')).parse_args()
  run(args)





