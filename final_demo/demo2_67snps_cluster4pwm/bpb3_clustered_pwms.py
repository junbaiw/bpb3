#this file is usded to evaluate bpb3 with clustered pwms by using bpb3 command lines
import os
import glob
from bpb3.script.script_high.other import common
import time

###########Setup main Path ################
#snp Data Input folder
source_data_folder= os.path.abspath('../../final_demo_data/demo1_67snps/in')
ref_seq = os.path.join(source_data_folder,"snps_ref.fasta")
alt_seq = os.path.join(source_data_folder,"snps_alt.fasta")

#Project Output path
project_name='demo2_67snps_cluster4pwm'

out_data_folder= os.path.abspath("../../final_demo_data/demo2_67snps_cluster4pwm")
common.check_folder(out_data_folder)

out_project_folder=os.path.join(out_data_folder, project_name)
common.check_folder(out_project_folder)

logging =common.setup_logger(out_project_folder,project_name, time.time())


# 1. Generate tempary pwms based on clustered and uncertain clustered PWMs from abc4pwm############
#here file folers 'abc4pwm_quality_assessed_out_seed5' and 'abc4pwm_uncertain_pwms' contain clustered pwms from abc4pwm by using 1772 pwms
#
in_clustered_pwms_path=out_data_folder
in_clustered_pwms_folder_name='abc4pwm_quality_assessed_out_seed5'
in_uncertain_pwms_folder_name='abc4pwm_uncertain_pwms'
out_tmp_pwms_folder_name='tmp_pwm'

#make pwm input files for clustred PWMs in bpb3 package
if True:
 logging.info("Make clustered PWMs input files for bpb3 ")
 cmd='bpb3 make_cluster4pwm --in_pwm_path ' + in_clustered_pwms_path +  \
        ' --in_clustered_pwm_folder ' +  in_clustered_pwms_folder_name + \
	' --out_file_folder ' + out_tmp_pwms_folder_name + \
	' --in_uncertain_pwm_folder ' + in_uncertain_pwms_folder_name 

 result_code=os.system(cmd)
 if result_code != 0:
      print("bpb3 make_cluster4pwm failed. See message in the console to find what went wrong.")
      exit(1)
 logging.info('Export at %s', os.path.join(in_clustered_pwms_path, out_tmp_pwms_folder_name))

######### Setup parameters in clustered PWMs for bayespi-bar3###################
# folder with PWM models for TF binding affinity, here we use clustered PWM first
pwm_folder=os.path.join(out_data_folder,out_tmp_pwms_folder_name)

# output folders  for bayespi-bar
foreground_folder = os.path.join(out_project_folder ,"foreground")
common.check_folder(foreground_folder)

block_res_folder = os.path.join(foreground_folder,"test_seq")
common.check_folder(block_res_folder)

# BayesPI-BAR options
# chemical potentials to use
chemical_potentials = "none -10 -13 -15 -18 -20"

# sequence shuffling iterations for dbA calculation
shuffling_iterations = 10000
random_seed=1

#Parallel option path
parallel_options_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "parallel_options.txt"))

# max rank of a PWM file in the foreground results to include in ranking TFs
max_ranking_for_TF = 50
top_ranked_TF_from_clusteredPWM=20
isChangePWMname=True
start_from_integration=False
export_top_ranked_TF=30
clean_tmp_files=False

#################START do bpb3 BayesPI-BAR on clustered PWMS ############
#2. run bpb3 on clustered PWMs - first level analysis
if True:
 logging.info("Pipeline parameters:")
 logging.info("\tBayesPI-BAR options: %d sequence shuffling iterations, %d chemical potentials: %s",
                 shuffling_iterations, len(chemical_potentials.split(" ")), chemical_potentials)

 cmd="bpb3 bayespi_bar" + \
      				" --reference_sequences " + ref_seq +\
                                " --alternate_sequences " + alt_seq +\
                                " --pwm_folder " + pwm_folder +\
                                " --chemical_potentials " + chemical_potentials +\
                                " --iterations " + str(shuffling_iterations) +\
                                " --result_folder " + block_res_folder +\
                                " --max_rank " + str(max_ranking_for_TF) +\
                                " --seed " + str(random_seed) + \
                                " --reuse_output" \
                                " @" + parallel_options_file
 result_code = os.system(cmd)
 if result_code != 0:
      logging.info("BayesPI-BAR failed. See message in the console to find what went wrong.")
      exit(1)
 logging.info('Export at %s', block_res_folder )

#collect all mlp in ranking results
#3. perform bayespi_bar on the selected pwms from the top ranked # TF from clusters of PWMs - second level analysis
output_file_folder='new_test_sep'
if True:
 cmd= 'bpb3 bpb3selectedPWM --cluster_out_path '  + out_data_folder + \
	' --bpb3_out_path ' + os.path.join(block_res_folder,'rankings')  + \
	' --in_fileString4clustered_pwm ' + in_clustered_pwms_folder_name  + \
	' --project_name ' + 'selected_pwm_from_cluster4rank_top' + str(top_ranked_TF_from_clusteredPWM) + \
	' --output_folder ' + output_file_folder  + \
	' --in_fileString4uncertain_pwm '+ in_uncertain_pwms_folder_name + \
	' --in_separateString4file '+ str("\'~\'") + \
	' --in_topRank_cutoff ' + str(top_ranked_TF_from_clusteredPWM) + \
	' --in_seq_ref_file ' + ref_seq +\
	' --in_seq_alt_file ' + alt_seq +\
	' --iterations ' + str(shuffling_iterations) +  \
	' --max_rank ' + str(export_top_ranked_TF) + \
	' --seed ' + str(random_seed) + \
	' --chemical_potentials ' + chemical_potentials +\
        ' --change_pwm_name ' * isChangePWMname + \
        ' --clean_tmp ' * clean_tmp_files + \
        ' --skip_dba_calculation ' * start_from_integration + \
	' --p_value_cutoff ' +'0.1' +  ' --integration pca' + \
         ' @' +  parallel_options_file
 result_code=os.system(cmd)
 if result_code !=0:
    logging.info("bpb3selectedPWM failed, see message in the console to find reason.")
    exit(1)
 logging.info('Export at %s ', out_data_folder )


if True:
 #6. clean temporary files
 logging.info("")
 logging.info("Remove temporary files in background or foreground calculations. ")
 foreground_folder1=os.path.join(out_data_folder, 'selected_pwm_from_cluster4rank_top' + str(top_ranked_TF_from_clusteredPWM), output_file_folder)
 background_folder= foreground_folder1
 background_pwm_folder=foreground_folder1
 print(foreground_folder1)
 result_code=os.system('bpb3 clean_tmp --background_out_file ' + background_folder + 
                ' --foreground_out_file ' + foreground_folder1  +  ' --chemical_potentials ' + chemical_potentials  +
                ' --background_pwm_file ' + background_pwm_folder + ' --clean_both 1' )
 if result_code !=0:
   logging.info("Error in clean tmp ")
   exit(1)

 foreground_folder1=foreground_folder
 print(foreground_folder1)
 result_code=os.system('bpb3 clean_tmp --background_out_file ' + background_folder +
                ' --foreground_out_file ' + foreground_folder1  +  ' --chemical_potentials ' + chemical_potentials  +
                ' --background_pwm_file ' + background_pwm_folder + ' --clean_both 1' )
 if result_code !=0:
   logging.info("Error in clean tmp ")
   exit(1)

else:
 logging.info('Skip step 6: remove temporary file %s ', background_pwm_folder) 

 


