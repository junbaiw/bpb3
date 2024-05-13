#this script is used to clean temporary files
import argparse
import sys
import glob 
import os
import pathlib

def my_parser(parser):
#parser= argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
#                                description=""" remove temporary files  """,
#                                add_help=False)

#try:
  required_named= parser.add_argument_group('required arguments')
  optional = parser.add_argument_group('optional arguments')

  required_named.add_argument('--background_out_file',help="Background calculation output file folder",
                              type=pathlib.Path, required=True, metavar='FILE FOLDER')
  required_named.add_argument('--foreground_out_file',help="Foreground calculation output file folder",
                              type=pathlib.Path, required=True, metavar="FILE FOLDER")
  required_named.add_argument('--background_pwm_file',help="Background PWM file folder",
                              type=pathlib.Path, required=True, metavar="FILE FOLDER")
  required_named.add_argument('--chemical_potentials', help='A list of chemical potentials used in the calculations, space separeted string',
                              type=str, required=True, metavar='Potentials',nargs="+")
  optional.add_argument('--clean_both',help="whether clean both foreground and background temporary files, default 0 for cleaning background, or 1 for foreground, or 2 for cleaning both foreground and background",
                        type=int, default=0, metavar="NUMBER")

#  optional.add_argument('-h','--help', action="help", help="show this help message and exit")
# 
#  args= parser.parse_args()
#except IOError as err:
#  print("Command line argument error: ", err, file=sys.stderr)
#  exit(1)
  return parser

def check_folder_file(out_file,file_prefer):
   """this function is used to check the exist of folder and return folder and its selected files in it"""
   f_file=os.path.abspath(out_file)
   if os.path.exists(f_file) and file_prefer != None :  
      all_f_file=glob.glob(os.path.join(f_file,file_prefer))
   else:
      all_f_file=None
   return f_file, all_f_file

def rm_file_and_folder(all_files,file_subfolders):
 if all_files !=None:
   for fi in all_files:
       #print(fi)
       if file_subfolders==None :
          cmd= 'rm -fr ' + fi
       elif "*" not in file_subfolders:
          cmd= 'rm -fr ' + os.path.join(fi,file_subfolders)
       else:
          tmp_files=os.path.join(fi,file_subfolders)
          cmd= 'rm -f ' + tmp_files
       print(cmd)
       result_code=os.system(cmd)
       if result_code !=0:
          print("Error in removing file folder ")
          exit(1)
 else:
   print(all_files, " not exist !")        

def run(args):
 foreground_file, all_foreground_files= check_folder_file(args.foreground_out_file,"*")
 background_file, all_background_files= check_folder_file(args.background_out_file,None)
 background_pwm_file, all_background_pwm_files= check_folder_file(args.background_pwm_file, "*")
 parallel_tmp_file, all_parallel_tmp_files= check_folder_file(args.foreground_out_file,"parallel_tmp")
 rm_file_and_folder(all_parallel_tmp_files, None)
     

 if args.clean_both ==0:
   print("Remove temporary files in background output: ")
   for chemical_potential in args.chemical_potentials:
     rm_file_and_folder([background_file], os.path.join("bayespi_bar_result","dba_alt","potential_"+chemical_potential ,"*0"))
     rm_file_and_folder([background_file], os.path.join("bayespi_bar_result","dba_ref","potential_"+chemical_potential ,"*0"))
   
   rm_file_and_folder([background_file], os.path.join("bayespi_bar_result","parallel_tmp"))
   rm_file_and_folder([background_pwm_file], None)
   print("\n")

 if args.clean_both ==1:
    print("Remove temporary files in foreround output: ")
    for chemical_potential in args.chemical_potentials:
       rm_file_and_folder(all_foreground_files,os.path.join("dba_alt","potential_"+chemical_potential,"*0"))
       rm_file_and_folder(all_foreground_files,os.path.join("dba_ref","potential_"+chemical_potential,"*0"))

    rm_file_and_folder(all_foreground_files,"parallel_tmp")

 if args.clean_both == 2:
   print("Remove temporary files in background output: ")
   #rm_file_and_folder([background_file],None)
   for chemical_potential in args.chemical_potentials:
     rm_file_and_folder([background_file], os.path.join("bayespi_bar_result","dba_alt","potential_"+chemical_potential,"*0"))
     rm_file_and_folder([background_file], os.path.join("bayespi_bar_result","dba_ref","potential_"+chemical_potential,"*0"))

   rm_file_and_folder([background_file], os.path.join("bayespi_bar_result","parallel_tmp"))
   rm_file_and_folder([background_pwm_file], None )
   print("\n")

   print("Remove temporary files in foreground output: ")
   for chemical_potential in args.chemical_potentials:
     rm_file_and_folder(all_foreground_files,os.path.join("dba_alt","potential_"+chemical_potential,"*0"))
     rm_file_and_folder(all_foreground_files,os.path.join("dba_ref","potential_"+chemical_potential,"*0"))

   rm_file_and_folder(all_foreground_files,"parallel_tmp")
   print("\n")
 print("Done with clean temporary files")

if __name__=='__main__':
  args=my_parser(argparse.ArgumentParser('python clean_tmp.py')).parse_args()
  run(args)
 














