#this file is usded to test bpb3 on 67 knonw regulation mutations that collected from
#2 Papers and HGMD database:
#https://pubmed.ncbi.nlm.nih.gov/19641089/
#https://pubmed.ncbi.nlm.nih.gov/25383969/
#https://pubmed.ncbi.nlm.nih.gov/24077912/
#which can be used to evauate the prediction accuracy of new bpb3
#
#some early applications of these 67 snps can be found in papers:
#https://pubmed.ncbi.nlm.nih.gov/26202972/
#https://pubmed.ncbi.nlm.nih.gov/31001324/
#
#In this demo, we do not run the whole pipeline but only use some of bayepi_bar command functions from bpb3 package
# 
import os
import glob
from bpb3.script.script_high.other import common
import time

###########Setup main Path ################
source_data_folder= os.path.abspath('../../final_demo_data/demo1_67snps/in')

# 1. Input ref and alt sequences of 67 SNPs
ref_seq = os.path.join(source_data_folder,"snps_ref.fasta")
alt_seq = os.path.join(source_data_folder,"snps_alt.fasta")

# 2. Project Output path
out_data_folder = os.path.abspath("../../final_demo_data/demo1_67snps/out")
common.check_folder(out_data_folder)
project_name='demo1_67snps'

# 3. setup log file
logging =common.setup_logger(out_data_folder,project_name, time.time())

# 4. folder with PWM models for TF binding affinity, here we use clustered PWM first
pwm_folder = os.path.join("../../final_demo_data/pwm")

# 5. output folders
foreground_folder = os.path.join(out_data_folder, "foreground")
common.check_folder(foreground_folder)
block_res_folder = os.path.join(foreground_folder,"test_seq")
common.check_folder(block_res_folder)

############bpb3 pipe line paramters ##################
# 6. BayesPI-BAR options
# chemical potentials to use
chemical_potentials = "none -10 -13 -15 -18 -20"

# sequence shuffling iterations for dbA calculation
shuffling_iterations = 10000

#Parallel option path
parallel_options_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "parallel_options.txt"))

# max rank of a PWM file in the foreground results to include in the background model
max_ranking_for_background = 50

logging.info("Pipeline parameters:")
logging.info("\tBayesPI-BAR options: %d sequence shuffling iterations, %d chemical potentials: %s",
                 shuffling_iterations, len(chemical_potentials.split(" ")), chemical_potentials)

result_code = os.system("bpb3 bayespi_bar" +
                                " --reference_sequences " + ref_seq +
                                " --alternate_sequences " + alt_seq +
                                " --pwm_folder " + pwm_folder +
                                " --chemical_potentials " + chemical_potentials +
                                " --iterations " + str(shuffling_iterations) +
                                " --result_folder " + block_res_folder +
                                " --max_rank " + str(max_ranking_for_background) +
                                " --reuse_output" 
                                " @" + parallel_options_file)

if result_code != 0:
      logging.info("BayesPI-BAR failed. See message in the console to find what went wrong.")
      exit(1)
logging.info('Export at %s', block_res_folder )


if True:
 #6. clean temporary files
 logging.info("")
 logging.info("Remove temporary files in background or foreground calculations. ")
 background_folder= foreground_folder
 background_pwm_folder=foreground_folder
 print(foreground_folder)
 result_code=os.system('bpb3 clean_tmp --background_out_file ' + background_folder + 
                ' --foreground_out_file ' + foreground_folder  +  ' --chemical_potentials ' + chemical_potentials  +
                ' --background_pwm_file ' + background_pwm_folder + ' --clean_both 1' )
 if result_code !=0:
   logging.info("Error in clean tmp ")
   exit(1)
else:
 logging.info('Skip step 6: remove temporary file %s ', background_pwm_folder) 





