#this script is usded to generate files for bayesPI-BAR3 based on DBD clustered pwms
#
import os
import glob
from bpb3.script.script_high.other import common
import argparse
#exec(open('test.py').read())

def my_parser(parser):
   required= parser.add_argument_group("Required ")
   required.add_argument("--in_pwm_path", help="File path contains all clustered PWMs and uncertain PWMs ",
                          required=True, type=str, metavar="FILE PATH")
   required.add_argument("--in_clustered_pwm_folder", help="File folder has all clustered PWMs, such as clustered_pwm_folder/DBD_name/out/cluster_number/repres/cluster_number_rep.mlp",
                          required=True, type=str, metavar="FOLDER")   

   optional= parser.add_argument_group("Option parameters , have default values")
   optional.add_argument("--out_file_folder", help="output file folder for renamed pwm files,default=tmp_pwm ",
                          type=str, metavar="FOLDER",default="tmp_pwm")
   optional.add_argument("--in_file_postfix", help="post prefix by file extension , default=*.mlp ", type=str, metavar="FILENAME EXTENSION", default= '*.mlp')
   optional.add_argument("--in_uncertain_pwm_folder", help="File folder has uncertain pwms or pwms are not assigned to clusters, default=uncertain_pwms",
                         type=str, metavar="FOLDER", default= 'uncertain_pwms')
   optional.add_argument("--in_separateString4file", help="string separation for joined new file name, default= ~", 
                       type=str,metavar="string separation", default='~')
   optional.add_argument("--in_folder_depth4clustered_pwm", help="folder depth for clustered pwms that will be used as new file name, default =-6",
                          type=int, metavar="FILE Depth", default=-6)
   optional.add_argument("--in_folder_depth4uncertain_pwm", help="folder depth for uncertain pwms that will be used as new file name, default = -2",
                           type=int, metavar="FILE Depth", default=-2)
   return parser


def copy_file2folder(in_folder1, files_in_folder1, out_folder,str4new_file,folder_depth):
  '''copy files to a tempory folder such as 'tmp_pwms' by using clustered PWMs. There are two input folder:
     in_folder1, a folder contains PWMs donot be assigned to any clusters after clustering quality check
     files_in_folder, all files in the folder1
     out_folder, is the output foler 
     str4new_file, new string to join file name such as '~'
    folder_depth, the last N names will be used as output file such as "uncertain_pwms" is -2 such as uncertain_pwms/12_SCRT1_1_from_SCRT1_jolma_DBD_M26-C2H2_ZF.mlp, 
     but DBD clustered PWMs is -6 such as 'quality_assessed_out_seed3/HMG/out/1/repres/1_rep.mlp' 
  '''
  for fi in files_in_folder1:
     fi=fi.replace('(','\(')
     fi=fi.replace(')','\)')
     out_file=str4new_file.join(fi.split(os.sep)[folder_depth:])
     #print(out_file)
     cmd = 'cp -f ' + os.path.abspath(fi) + ' ' + os.path.join(os.path.abspath(out_folder) , out_file)
     result_code=os.system(cmd)
     if result_code !=0:
        print('Error in ', cmd)
        exit(1)
  print('Copied ', str(len(files_in_folder1)),' from ', in_folder1, ' to ', os.path.abspath(out_folder) )

def run(args):
  in_pwm_folder=args.in_pwm_path
  in_folder1=os.path.join(in_pwm_folder,args.in_clustered_pwm_folder)
  in_folder2=os.path.join(in_pwm_folder,args.in_uncertain_pwm_folder)
  in_file_postfix=args.in_file_postfix
  out_folder=os.path.join(in_pwm_folder,args.out_file_folder)
  #common.check_folder(out_folder)
  common.prepare_result_folder(out_folder)

  #get file names, here use default output file structure from abcpwm
  #here clustered file structure is fixed such as quality_assessed_out_seed3/HMG/out/1/repres/1_rep.mlp
  files_in_folder1=glob.glob(os.path.join(in_folder1,'*','out','*','repres',in_file_postfix))
  files_in_folder2=glob.glob(os.path.join(in_folder2,in_file_postfix))

  #copy files
  copy_file2folder(in_folder1,files_in_folder1, out_folder, 
                   str4new_file=args.in_separateString4file, folder_depth=args.in_folder_depth4clustered_pwm)
  copy_file2folder(in_folder2, files_in_folder2, out_folder,
                  str4new_file=args.in_separateString4file, folder_depth=args.in_folder_depth4uncertain_pwm)
  return


if __name__=='__main__':
  args= my_parser(argparse.ArgumentParser('python make_cluster4pwm.py')).parse_args()
  run(args)


