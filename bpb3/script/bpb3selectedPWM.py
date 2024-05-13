#this script is used to run bayespi-bar for snp with selected pwms which
#means each snp/block has difffedent input pwms . For example, select the top ranked pwm/cluster from the first round application of 
#bpb3 on DBD clustered PWMs
import glob
import pandas as pd
import os
from bpb3.script.script_high.other import common
from bpb3.script.script_high import parallel, bayespi_common
from bpb3.script.script_high.other import check_accuracy as ck_ac
from bpb3.script.bayespi_bar import my_parser, collect_dba, calculate_delta_dba, integrate_delta_dba, compute_rankings
import argparse
import multiprocessing as mp
import math
import shutil

#exec(open("bayespi_bar4selectedPWM.py").read())
def my_parser2(parser):
  #here is new parameters for clustered PWM in bayespi-bar calculation
   required = parser.add_argument_group('required arguments for clustering')

   required.add_argument('--cluster_out_path', help="File path for clustered PWMs results", required=True,
                              type=str, metavar="FILE PATH")
   required.add_argument('--bpb3_out_path', help="File path of bpb3 ranking based on clustered PWMs", required=True,
                              type=str, metavar="FILE PATH")
   required.add_argument("--in_fileString4clustered_pwm", help="File string or file folder/name used/contained for clustered PWMs such as quality_assessed_out_seed3",
                          required=True, type=str, metavar="FILE String")
   
   required.add_argument('--project_name',help="project name for the current calculation", required=True,
                            type=str, metavar="PROJECT NAME")
   required.add_argument('--output_folder', help="Output folder name for the current calculation", required=True,
                           type=str, metavar="FILE FOLDER NAME")
   required.add_argument('--in_seq_ref_file', help="File name of reference genome sequence file", required=True,
                            type=str, metavar="FILE NAME")
   required.add_argument('--in_seq_alt_file', help="File name of alternative genome sequence file", required=True,
                            type=str, metavar="FILE NAME")


   optional = parser.add_argument_group('optional arguments for clustering')
   optional.add_argument("--in_fileString4uncertain_pwm", help="File string or file folder/name used/contained for untertain PWMs, default is uncertain_pwms",
                        type=str, metavar="STRING NAME", default='uncertain_pwms')
   optional.add_argument("--in_separateString4file", help="string separation for joined new file name, default= \'~\'", 
                       type=str,metavar="string separation", default='~')
   optional.add_argument("--in_topRank_cutoff", help="The cutoff value of top ranked results from the first level analysis of clustered PWMs , default=10 ",
                       type=int, metavar="NUMBER", default=10)
   optional.add_argument("--skip_dba_calculation", help="Skip dba calculation in bbp3 if the results are already exist, default is False",
                         action="store_true")
   optional.add_argument("--clean_tmp", help ="Clean temporary files in snp calculations default = false ", action="store_true")   
   optional.add_argument("--export_selected_pwms", help="Export all selected pwms from cluster4pwm to a folder: selected_pwm default =false", action="store_true")
   optional.add_argument("--export_pwm_to_different_folder", help="Export folder path for pwms that are selected from cluster4pwm , default=false and all pwms will be exported to the folder: selected_pwm ",
                          action="store_true")
   optional.add_argument("--change_pwm_name", help="Change pwm names for results of cluster4pwm by removing the first few labels of clustering information in each PWM, default =False ", action="store_true")
   #add bayespi_bar parser
   my_parser(parser)

   return parser

def find_ddba4blocks(out_data_folder, potentials):
 #for block_name in significant_block_tags:
   #all_snp_files=glob.glob(os.path.join(in_cluster4pwm_folder,block_name,'snp_*'))
   all_snp_files=glob.glob(os.path.join(out_data_folder, 'snp_*'))
   #print(len(all_snp_files))
   all_snp_files2=[ i for i in all_snp_files if os.path.basename(i).split('_')[1].isnumeric()]
   #print(all_snp_files2)
   #print(len(all_snp_files2))
   record_potentials={}
   record_ddba_mean=[]
   record_ddba_pca=[]

   for snp in all_snp_files2:
      tmp_folder=os.path.join(snp,'ddba')
      for pot in potentials:
          tmp_potential='potential_'+pot+'_ddba.tsv'
          if tmp_potential in record_potentials.keys():
             record_potentials[tmp_potential].append(pd.read_csv(os.path.join(tmp_folder,tmp_potential),sep='\t').copy())
          else:
             record_potentials[tmp_potential]=[pd.read_csv(os.path.join(tmp_folder,tmp_potential),sep='\t').copy()]
      record_ddba_mean.append(pd.read_csv(os.path.join(tmp_folder,'ddba_integrated_using_mean_ddba.tsv'),sep='\t').copy())
      record_ddba_pca.append(pd.read_csv(os.path.join(tmp_folder,'ddba_integrated_using_pca.tsv'), sep='\t').copy())
   return record_potentials, record_ddba_mean, record_ddba_pca, all_snp_files2
   

def merge_dfs(t_ddba):
  tm=pd.DataFrame()
  for ti in t_ddba:
    if tm.empty:
      tm=ti.copy()
    else:
      tm=tm.merge(ti,how='outer',suffixes=('_x','_y'),left_on='Sequence_name',right_on='Sequence_name' ).copy()
  tm=tm.fillna(0)
  return tm.copy()

def concat_dfs(t_ddba):
  tm=pd.DataFrame()
  for ti in t_ddba:
     ti=ti.set_index('Sequence_name')
     if tm.empty:
        tm=ti.copy()
     else:
        tm=pd.concat([tm,ti],axis=1).copy()
  tm=tm.fillna(0)
  return tm.copy()

def change_cluster4pwm_name(pwm_str):
 if 'Sequence_name' not in pwm_str:
   pwm1=pwm_str.split('_')
   #remove the first cluster number from PWM name
   if pwm1[0].isnumeric():
     pwm1='_'.join(pwm1[1:])
   else:
     pwm1=pwm_str
   pwm2='-'.join(pwm1.split('-')[0:-1])+'.mlp'
 else:
   pwm2=pwm_str
 return pwm2

def export_ddba2file(ddba_mean, out_file, isChangeName):
  col2mean=ddba_mean.columns.tolist()
  if isChangeName:
    new_col2mean=[ change_cluster4pwm_name(i) for i in col2mean ]
  else:
    new_col2mean=col2mean

  ddba_mean.columns=new_col2mean

  uq_ddba_mean=ddba_mean.groupby(level=0,axis=1).sum().copy()
  uq_ddba_mean=uq_ddba_mean.reset_index().copy()

  #out_file='ddba_mean.tsv'
  uq_ddba_mean.to_csv(out_file, sep='\t',index=False)
  return uq_ddba_mean


def list_of_list2uniqList(a_list_of_lists):
  uq_list=list(set([a for b in a_list_of_lists for a in b]))
  return uq_list

def find_filename4ranked_pwm(record_snp2rankings, quality_string, uncertain_string,str4file_sep,cluster_out_folder):
  #for each record_snp2rankings find its mlp file names that ranked in top N
  record_snp2ranked_mlps={}
  for ri in record_snp2rankings.keys():
    tmp_file_name = record_snp2rankings[ri].file_name.to_list()
    record_mlps=[]
    for ii in tmp_file_name:
      if ii.startswith(quality_string):
        tmp_path=os.path.split(os.path.split(os.path.join(os.path.abspath(cluster_out_folder), ii.replace(str4file_sep,'/')))[0])[0]
        tmp_mlp=glob.glob(os.path.join(tmp_path,'*.mlp'))
      elif ii.startswith(uncertain_string):
        tmp_mlp=os.path.join(os.path.abspath(cluster_out_folder),ii.replace(str4file_sep,'/'))
      else:
        print('Unknonw type cluster : ' , ii )
        pass
      if not isinstance(tmp_mlp,list):
          tmp_mlp=[tmp_mlp]
      record_mlps.append(tmp_mlp)
    record_snp2rankings[ri]['pwm_files']=record_mlps
    record_snp2ranked_mlps[ri]= list_of_list2uniqList(record_mlps )
  return record_snp2ranked_mlps

def make_seq_files(record_snp2ranked_mlps,out_data_folder,ref_seq,alt_seq):
  #for each snp extract its ref and alt sequences as well as create export folders
  #such as dba_ref, dba_alt, ddba and rankings, seqs
  record_id={}
  record_ref_seq={}
  record_alt_seq={}
  record_dba_folder={}
  block_id =0
  #export snp sequenes
  for ri in record_snp2ranked_mlps.keys():
    block_name='snp_'+ str(block_id)
    record_id[ri]=block_name
    tmp_path=os.path.join(out_data_folder,block_name)
    if not os.path.exists(tmp_path):
      os.mkdir(tmp_path)
    tmp_seq_folder=os.path.join(tmp_path,'seqs')
    if not os.path.exists(tmp_seq_folder):
      os.mkdir(tmp_seq_folder)

    ref_dba_folder=os.path.join(tmp_path,'dba_ref')
    alt_dba_folder=os.path.join(tmp_path,'dba_alt')
    ddba_folder=os.path.join(tmp_path,'ddba')
    rank_folder=os.path.join(tmp_path,'rankings')
    if not os.path.exists(ref_dba_folder):
       os.mkdir(ref_dba_folder)
    if not os.path.exists(alt_dba_folder):
       os.mkdir(alt_dba_folder)
    #common.prepare_result_folder(ref_dba_folder)
    #common.prepare_result_folder(alt_dba_folder)
    #common.prepare_result_folder(ddba_folder)
    #common.prepare_result_folder(rank_folder)
    #record_dba_folder contains, ref_dba, alt_dba, ddba, and rank folders for each snp 
    record_dba_folder[ri]=(ref_dba_folder, alt_dba_folder, ddba_folder,rank_folder)

    tmp_ref_seq= [ii for ii in ref_seq if ii[0]==ri]
    tmp_alt_seq= [ii for ii in alt_seq if ii[0]==ri]
    tmp_ref_file=os.path.join(tmp_seq_folder,block_name+'_ref.fasta')
    tmp_alt_file=os.path.join(tmp_seq_folder,block_name+'_alt.fasta')
    common.write_fasta(tmp_ref_seq,tmp_ref_file)
    common.write_fasta(tmp_alt_seq,tmp_alt_file)

    record_ref_seq[ri]=tmp_ref_file
    record_alt_seq[ri]=tmp_alt_file
    block_id +=1
  return record_ref_seq, record_alt_seq, record_id, record_dba_folder

def make_dba_commands(chemical_potentials, res_folder,record_snp2ranked_mlps, record_seq,record_id, iterations, seed, str4subfolder):
 #make bayespi-bar commmad line for bpb3 ref_seq or alt_seq calculations
 commands=[]
 for potential in chemical_potentials:
   #potential_result_folder=os.path.join(res_folder, "potential_" + potential)
   for ri in record_snp2ranked_mlps.keys():
      tmp_pwm_files=record_snp2ranked_mlps[ri]
      tmp_seq_file=record_seq[ri]
      tmp_folder=os.path.join(record_id[ri],str4subfolder)
      potential_result_folder=os.path.join(res_folder, tmp_folder,'potential_'+potential)
      common.prepare_result_folder(potential_result_folder)
      for pwm_fi in tmp_pwm_files:
         commands.append(bayespi_common.compute_affinity_command(tmp_seq_file, pwm_fi, iterations,
                                                                    potential_result_folder, potential, seed))
 return commands


def collect_all_ranked_results(in_file_folder, record_id, isChangeName):
  #collect ranked results from all snps to one folder
  #in_file_folder='../../data/cluster4pwm_old/ranked_pwms4snp/new_test/'
  out_file_folder=os.path.join(in_file_folder,'rankings')
  if not os.path.exists(out_file_folder):
     print('Create folder: ', out_file_folder)
     os.mkdir(out_file_folder)
     #os.system('ln -s '+ out_file_folder + ' '+ 'rankings')

  snp_ids=[]
  for ke, val in record_id.items() :
    snp_id=val
    in_ranking_file=glob.glob(os.path.join(in_file_folder,snp_id,'rankings','*.tsv'))
    for fi in in_ranking_file:
      #cmd='cp -f ' + fi + ' ' + out_file_folder
      #os.system(cmd)
      #remove cluster number and DBD before copy file to a new folder
      fi_df=pd.read_csv(fi,sep='\t',skiprows=1)
      #read the first line of comment in file
      first_line=pd.read_csv(fi,nrows=0)
      if isChangeName:
        fi_df["file_name"]=fi_df.file_name.apply(lambda x: change_cluster4pwm_name(x))
      out_file=os.path.join(out_file_folder, os.path.basename(fi) )
      #write first line comment to file
      with open(out_file,'w') as of:
          of.write(first_line.columns.to_list()[0]+'\n')
      of.close()
      fi_df.to_csv(out_file,sep='\t',index=False,mode='a')
      #print(out_file)
  print('Export all ranked results at: ', out_file_folder)

def collect_all_ddba_results(in_data_folder, potentials, isChangeName ):
   out_file_folder=os.path.join(in_data_folder,'ddba')
   if not os.path.exists(out_file_folder):
     print('Create folder: ', out_file_folder)
     os.mkdir(out_file_folder)

   t_potential,t_ddba_mean, t_ddba_pca, t_all_snp_files=find_ddba4blocks(in_data_folder, potentials)
   #for bk_key in record_blocks.keys():
   #  t_potential,t_ddba_mean, t_ddba_pca=record_blocks[bk_key]
   ddba_mean=concat_dfs(t_ddba_mean)
   ddba_pca=concat_dfs(t_ddba_pca)

   uq_ddba_mean=export_ddba2file(ddba_mean,os.path.join(out_file_folder,'ddba_integrated_using_mean_ddba.tsv'), isChangeName)
   uq_ddba_pca=export_ddba2file(ddba_pca,os.path.join(out_file_folder,'ddba_integrated_using_pca.tsv'), isChangeName)
     
   #print(uq_ddba_mean.shape)
   #print(uq_ddba_pca.shape)

   for pt in t_potential.keys():
       tmp_df=t_potential[pt]
       pt_df=concat_dfs(tmp_df)
       uq_pt=export_ddba2file(pt_df, os.path.join(out_file_folder,pt), isChangeName)
       #print(uq_pt.shape)
   return t_all_snp_files

def collect_calculate_integrate_dba4ranking( args):
   record_dba_folder, block_chunks, targs= args
   for ki in block_chunks:
      collect_dba(record_dba_folder[ki][0])
      collect_dba(record_dba_folder[ki][1])
      calculate_delta_dba(record_dba_folder[ki][0], record_dba_folder[ki][1], record_dba_folder[ki][2], targs.p_value_cutoff, targs.normalize_dba)
      integrate_delta_dba(record_dba_folder[ki][2], record_dba_folder[ki][2], "mean_ddba", targs.chemical_potentials)
      integrate_delta_dba(record_dba_folder[ki][2],record_dba_folder[ki][2],targs.integration, targs.chemical_potentials)
      compute_rankings(record_dba_folder[ki][2], record_dba_folder[ki][3],targs.integration, targs.max_rank, targs.med_ddba_p_value)

def run(args):
  #read all bpb3 results 
  in_files=glob.glob(os.path.join(os.path.abspath(args.bpb3_out_path),'*.tsv'))
  if len(in_files)<1:
     print('Input bpb3 ranked snps not find, I STOP! Please check folder, ', args.bpb3_out_path)
     exist(1)

  record_snp2rankings=ck_ac.load_results(in_files, rankCutoff=args.in_topRank_cutoff)
  print('Select top '  , args.in_topRank_cutoff, ' ranked TFs for calculation!')
  record_snp2ranked_mlps=find_filename4ranked_pwm(record_snp2rankings, args.in_fileString4clustered_pwm, 
                   args.in_fileString4uncertain_pwm,args.in_separateString4file,args.cluster_out_path)

  out_data_folder0=os.path.join(os.path.abspath(args.cluster_out_path),args.project_name)
  out_data_folder=os.path.join(out_data_folder0, args.output_folder)
  if not args.skip_dba_calculation :
     if not os.path.exists(out_data_folder0):
         print('Create project folder: ', out_data_folder0)
         os.mkdir(out_data_folder0)
     else:
          print('Project folder exists', out_data_folder0)
          #common.clear_folder(out_data_folder0)

     if not os.path.exists(out_data_folder):
         print('Create output folder: ', out_data_folder)
         os.mkdir(out_data_folder)
     else:
         print('Output folder exists, clean it now' , out_data_folder)
         common.clear_folder(out_data_folder)
  else:
     print('Skip ddba calculation!')

  ref_seq=common.read_fasta(os.path.abspath(args.in_seq_ref_file))
  alt_seq=common.read_fasta(os.path.abspath(args.in_seq_alt_file))

  #make seq files for bayespi-bar calculation
  record_ref_seq, record_alt_seq, record_id, record_dba_folder= make_seq_files(record_snp2ranked_mlps,out_data_folder,ref_seq,alt_seq)
  commands=[] 
  if not args.skip_dba_calculation :
    #make commmad line for bpb3 ref_seq and alt_seq
    commands.extend(make_dba_commands(args.chemical_potentials, out_data_folder,
                                             record_snp2ranked_mlps, record_ref_seq, record_id,args.iterations, args.seed,'dba_ref'))
    commands.extend(make_dba_commands(args.chemical_potentials, out_data_folder,
                                             record_snp2ranked_mlps, record_alt_seq, record_id,args.iterations, args.seed,'dba_alt'))
 
    runner = parallel.Parallel(args, os.path.join(out_data_folder, "parallel_tmp"))
    runner.run_commands(commands,"dba")
  
  #equality distribute record_id keys to num_of_process for computing ddba and ranking
  num_of_processes=args.max_nodes
  all_records_dba=list(record_dba_folder.keys())
  num_in_chunks=int(math.ceil(len(all_records_dba)/num_of_processes))
  block_chunks=[ all_records_dba[x:x+num_in_chunks]  for x in range(0, len(all_records_dba),num_in_chunks)]
  print(block_chunks[-1])
  #print(record_dba_folder.keys())
  #print(len(block_chunks))

  if num_of_processes > len(block_chunks):
     num_of_processes= len(block_chunks)
  pool=mp.Pool(processes= num_of_processes )
  pool.map(collect_calculate_integrate_dba4ranking, [(record_dba_folder,block_chunks[loop],args) for loop in range(0, num_of_processes)],1)

  pool.close()
  pool.join()

  #export all results in one folder
  collect_all_ranked_results(out_data_folder, record_id, args.change_pwm_name)
  all_snp_files=collect_all_ddba_results(out_data_folder,  args.chemical_potentials, args.change_pwm_name )

  #export record_id for snps
  #Export snp id to snp file
  out_df=pd.DataFrame(list(record_id.items()))
  out_file=os.path.join(out_data_folder,'snp_id2snp.tsv')
  out_df.to_csv(out_file, sep='\t',index=False, header=None)

  if args.export_selected_pwms:
     #print(record_snp2ranked_mlps)
     #uq_record_snp2ranked_mlps=list(set(record_snp2ranked_mlps))
     if not args.export_pwm_to_different_folder:
       #export selected pwms of cluster4pwm to a common folder selected_pwms
       out_pwm_folder=os.path.join(args.cluster_out_path,"selected_pwm")
       if not os.path.exists(out_pwm_folder):
         print("Create folder for selected PWMs from cluster4pwm,  "  , out_pwm_folder)
         os.mkdir(out_pwm_folder)

       all_recorded_mlps=[]
       for bk_key in record_snp2ranked_mlps.keys():
         all_recorded_mlps=all_recorded_mlps+ list(set(record_snp2ranked_mlps[bk_key]))
       uq_record_snp2ranked_mlps=list(set(all_recorded_mlps))
       for pwi in uq_record_snp2ranked_mlps:
            if ')' in pwi or '(' in pwi:
               pwi=pwi.replace(')','\)')
               pwi=pwi.replace('(','\(')

            os.system('\cp -fr ' + pwi + ' ' +  out_pwm_folder)
       print("copied " +  str(len(uq_record_snp2ranked_mlps)) + ' pwm files')
     else:
       for bk_key in record_snp2ranked_mlps.keys():
           #bk_id=bk_key.split('_patient_')[0]
           out_pwm_folder=os.path.join(out_data_folder, "selected_pwm")
           if not os.path.exists(out_pwm_folder):
              print("Create folder for selected PWMs from cluster4pwm, ", out_pwm_folder)
              os.mkdir(out_pwm_folder)
           for pwi in record_snp2ranked_mlps[bk_key]:
              if ')' in pwi or '(' in pwi:
                pwi=pwi.replace(')','\)')
                pwi=pwi.replace('(','\(')
              os.system('\cp -fr ' + pwi + ' ' +  out_pwm_folder)
           print("copied " +  str(len(record_snp2ranked_mlps[bk_key])) + ' pwm files')


  #remove tempartory files
  if args.clean_tmp:
    print('Removing temporary files in cluster4pwm : ' + out_data_folder)
    for si in all_snp_files:
      shutil.rmtree( os.path.join(out_data_folder, si) ) 
    shutil.rmtree(os.path.join(out_data_folder, 'parallel_tmp')) 

  return commands

if __name__== "__main__":
##################
 #load args from bpb3
 args=my_parser2(argparse.ArgumentParser('python bpb3selectedPWM.py')).parse_args()
 commands=run(args)


 if False:
  #args for pwm clustering
  #args.cluster_out_folder='../../data/cluster4pwm_sep/'
  args.cluster_out_folder='../data/out/sep7_run2/'

  #bpb3 for PWM clustering results
  #args.bpb3_out_folder='../../data/test_67snps/out_sep/foreground/test_seq/rankings/'
  args.bpb3_out_folder='../data/out/sep7_run2/test_67snps/foreground/test_seq/rankings/'

  args.quality_string='quality_assessed_out_seed3'
  args.uncertain_string='uncertain_pwms'
  args.str4file_sep='~'
  #exported path for ranks and bpb3 results
  args.project_name='selected_pwm_from_cluster4rank_top20'
  args.output_folder='new_test_sep'

  #selection of top ranked TF from cluster pwm results
  args.topRankedTFs=20

  #skip dba calculations for ref_seq and alt_seq
  args.is_skip=False

  #below are args in bpb3
  args.in_seq_ref_file='../../data/test_67snps/in/snps_ref.fasta'
  args.in_seq_alt_file='../../data/test_67snps/in/snps_alt.fasta'
  
  args.use_cores=10
  args.seed=1
  args.iterations=10000
  args.chemical_potentials=['-8','-10','-13','-15','-18','-20']
  args.use_slurm=True
  args.slurm_account='nn4605k'
  args.max_nodes=10
  args.p_value_cutoff=0.1
  args.normalize_dba= False
  args.integration='pca'
  args.max_rank=30
 
  commands=run(args)







