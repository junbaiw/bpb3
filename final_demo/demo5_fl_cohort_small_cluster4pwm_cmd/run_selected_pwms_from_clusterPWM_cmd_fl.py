import os
from bpb3.script.script_high.other import common
import time


#This script is used to evaluate effect of TF binding affinity changes due to DNA sequence mutations based on 
#selected PWMs from clustered PWMs by abc4pwm - second level analysis
#Here, we assume all mutation blocks and the first level analysis based on clustered PWMs/representative motifs from abc4pwm  are already done.
#This script simply selects a set of top ranked TFs/clustered PWMs from the first level analysis, then use all PWMs clustered in these seected clusters for
#the second level analysis, to find significant TF binding affinity changes. 
#Here, most of bpb3 command line functions and parameter configures are shown in the same file
#

#parameter from bpb3
in_data_folder='../../final_demo_data/demo5_fl_cohort_small_cluster4pwm_cmd/'
mussd_out_data_folder=os.path.join(in_data_folder, 'out')
significant_blocks_file = os.path.join(mussd_out_data_folder, "significant_blocks.txt")
shuffling_iterations=10000
random_seed=1
chemical_potentials= "none -10 -13 -15 -18 -20"
parallel_options_file= os.path.abspath(os.path.join(os.path.dirname(__file__), "parallel_options.txt"))
significant_block_tags = frozenset(common.read_lines(significant_blocks_file))

#Parameter for cluster4pwm
in_clustered_pwms_folder_name='abc4pwm_quality_assessed_out_seed5'
in_uncertain_pwms_folder_name='abc4pwm_uncertain_pwms'
#select top ranked TFs from clusters this one seems not work?? 
top_ranked_TF_from_clusteredPWM=25

#project name or file folder for storing calculated results.
project_name='bpb3_by_selected_pwm_from_cluster4pwm'
out_project_folder=os.path.join(in_data_folder, project_name)

#parameter for tumor or normal expression data
group1_str='tumor'
group2_str='normal'
#Export file path
group1_median_rpkm_file = os.path.join(mussd_out_data_folder, group1_str + "_median_patient_rpkm.tsv")
group2_median_rpkm_file = os.path.join(mussd_out_data_folder, group2_str +"_median_patient_rpkm.tsv")
#parameter for bpb3 mussd
# output folders for mussd
mussd_result_folder=os.path.join(mussd_out_data_folder, 'mussd_blocks')

#parameter for new selected PWMs from clusters
foreground_folder = os.path.join(in_data_folder,project_name )
background_folder = os.path.join(in_data_folder, "background")
background_pwm_folder = os.path.join(in_data_folder, "pwm_for_background")

#parameter for genpme information and background 
background_command_file= os.path.join(in_data_folder, "make_background.sh")
genome_folder=os.path.abspath("../../final_demo_data/genome_data/hg19")
genome_fasta_file =os.path.join(genome_folder, "human_g1k_v37_decoy.fasta") 
regions_file= os.path.join(mussd_out_data_folder, "regions.bed")

#path for pwms selected from all clusters , this selected_pwm folder is hard coded in script
selected_pwm_folder=os.path.join(in_data_folder,'selected_pwm')

#background parameters
background_shuffling_iterations=1000
resampling_iterations=10
gene_samples=8
mutation_signature_file=None
tumor_mutations_as_background=True
max_ranking_for_background=15

#signifcant block parameters
max_ranking_for_significant_blocks= 25
affinity_change_p_value=0.01
pval_correction4affinity=False


logging=common.setup_logger(in_data_folder, project_name, time.time())

##assume both forground and background ddba are already exists!
#this  works for foreground as skip_dba_calculation  
start_from_integration=False

runBayesPiBar4selectedPwms=True
#start from intergration only works here!
runBackgroundModel=True
runBlockSignificance=True
filterTFbyGene=True
runRemoveTempFile=True

############Below are calculations for PWMs selected from clustered PWMs###############
if runBayesPiBar4selectedPwms:
 #loop in each block for computing ddba
 for block_name in significant_block_tags :
  ref_seq=os.path.join(mussd_out_data_folder,'mussd_blocks/'+block_name + '_ref.fasta')
  alt_seq=os.path.join(mussd_out_data_folder,'mussd_blocks/'+block_name + '_alt.fasta')
  block_res_folder=os.path.join(mussd_out_data_folder,'foreground/'+ block_name ,'rankings')

  #collect all mlp in ranking results
  #1. perform bayespi_bar on selected pwms from top ranked TF clusters
  logging.info("")
  logging.info("Step 1: computing the foreground model for selected PWMs from clusters")
  cmd= 'bpb3 bpb3selectedPWM --cluster_out_path '  + in_data_folder + \
        ' --bpb3_out_path ' + block_res_folder  + \
        ' --in_fileString4clustered_pwm ' + in_clustered_pwms_folder_name  + \
        ' --project_name ' + project_name + \
        ' --output_folder ' + block_name + \
        ' --in_fileString4uncertain_pwm '+ in_uncertain_pwms_folder_name + \
        ' --in_separateString4file '+ str("\'~\'") + \
        ' --in_topRank_cutoff ' + str(top_ranked_TF_from_clusteredPWM) + \
        ' --in_seq_ref_file ' + ref_seq +\
        ' --in_seq_alt_file ' + alt_seq +\
        ' --iterations ' + str(shuffling_iterations) +  \
        ' --max_rank ' + str(top_ranked_TF_from_clusteredPWM) + \
        ' --seed ' + str(random_seed) + \
        ' --clean_tmp ' * runRemoveTempFile + \
        ' --skip_dba_calculation ' * start_from_integration + \
        ' --chemical_potentials ' + chemical_potentials +\
        ' --export_selected_pwms ' + \
        ' --p_value_cutoff ' +'0.1' +  ' --integration pca' + \
         ' @' +  parallel_options_file
  result_code=os.system(cmd)
  if result_code !=0:
    logging.info("bpb3selectedPWM failed, see message in the console to find reason.")
    exit(1)
  logging.info('Export at %s , %s ', os.path.join(in_data_folder,project_name,block_name), block_name )
else:
  logging.info("Skip step 1: computing the foreground model for selected PWMs from clusters")

if runBackgroundModel:
   # 2. Compute the background for PWMs selected from clustered PWMs
   logging.info("")
   logging.info("Step 2: Computing the background model for selected PWMs from clusters")
   logging.info("Creating shell file with the command to compute the background: %s", background_command_file)
   result_code = os.system("bpb3 choose_background_parameters " +
                            " --blocks_folder " + mussd_result_folder +
                            " --foreground_folder " + foreground_folder +
                            " --background_pwm_folder " + background_pwm_folder +
                            " --output_file " + background_command_file +
                            " --block_ids " + " ".join(significant_block_tags) +
                            " --genome " + genome_fasta_file +
                            " --regions " + regions_file +
                            " --pwm_folder " + selected_pwm_folder +
                            " --iterations " + str(background_shuffling_iterations) +
                            " --chemical_potentials " + chemical_potentials +
                            " --background_output_folder " + background_folder +
                            " --block_resample_iterations " + str(resampling_iterations) +
                            " --block_samples_to_take " + str(gene_samples) +
                            " --reuse_output" + 
                            " --start_from_integration" * start_from_integration +  
                            " --max_rank " + str(max_ranking_for_background) +
                            ((" --mutation_signature " + mutation_signature_file)
                             if mutation_signature_file is not None else "") +
                            (" --background_mutations %s" %
                             os.path.join(mussd_result_folder, "mutations_in_regions.tsv")) *
                            tumor_mutations_as_background +
                            " @" + parallel_options_file)

   if result_code !=0:
       logging.info("Error in choose background parameters ")
       exit(1)

   logging.info('Export at %s ', background_command_file )

   if os.path.exists(background_folder):
      logging.info("The folder with background output already exists. The background computation script will now "
                     "check its consistency and recompute the broken or missing data")
   else:
      logging.info("The background computation script will now compute the background delta-dbA score distributions %s ", background_command_file)

   result_code = os.system("bash " + background_command_file)
   if result_code != 0:
     logging.info("Background model computation failed. See messages in the console to find what went wrong.")
     exit(1)
 
   logging.info("Export at %s ", background_folder)

else: 
 logging.info('Skip step 2 background_computation %s ', background_folder)



if runBlockSignificance:
 # 3. Run significance tests for TF affinity changes between foreground and background model
 logging.info("")
 logging.info("Step 3: Significance tests will now be computed comparing background and foreground delta-dbA "
                 "distributions for each PWM in each mutation block")

 for block_tag in significant_block_tags:
        block_folder = os.path.join(foreground_folder, block_tag)
        result_file = os.path.join(block_folder, "result.tsv")
        print( block_tag, ":" )
        result_code = os.system("bpb3 affinity_change_significance_test " +
                                " --background_folder " + os.path.join(background_folder, "bayespi_bar_result") +
                                " --foreground_folder " + block_folder +
                                " --output_file " + result_file +
                                " --max_rank " + str(max_ranking_for_significant_blocks) +
                                " --p_value " + str(affinity_change_p_value) +
                                " --pval_correction " * pval_correction4affinity)
        if result_code !=0:
            logging.info("Error in affinity change significance test %s", block_tag)
            exit(1)
        
        logging.info("The table of significantly affected PWMs for block %s is stored in the file %s ", block_tag,
                     result_file)
else:
  logging.info("Skip step 3: significance tests for comparing background and foreground delta.dba distribution for each PWM in each mutation block")

if filterTFbyGene:
 # 4. Filter results by TF gene expression data if there is one
 logging.info("")
 logging.info("Step 4: PWMs corresponding to the non-expressed TFs will now be removed from the results")
 for block_tag in significant_block_tags:
        block_folder = os.path.join(foreground_folder, block_tag)
        result_file = os.path.join(block_folder, "result.tsv")
        filtered_result_file = os.path.join(block_folder, "result_filtered.tsv")

        print("Filtering results by gene expression for block", block_tag, "...")
        result_code = os.system("bpb3 filter_results_by_gene_expression_cluster4pwm " +
                                " --results_file " + result_file +
                                " --expression " + group1_median_rpkm_file +
                                " --notSepString4cluster " + 
                                " --output_file " + filtered_result_file)

        if result_code != 0:
            logging.info("Filtering results for block %s failed. "
                         "See messages in the console to find what went wrong.", block_tag)
            exit(1)

        logging.info("The table of significantly affected PWMs for block %s, containing only those TFs which are "
                     "expressed in patients, is stored in the file %s", block_tag,
                     filtered_result_file)

 #5. plot all heatmaps
 logging.info("")
 logging.info("Step 5: Generate plots or heatmaps for the results")

 result_code=os.system('bpb3 make_plots --foreground_folder ' + foreground_folder +
          ' --background_folder ' + background_folder + 
          ' --mussd_result_folder ' + mussd_result_folder )

 if result_code !=0:
       logging.info("Error in make plots ")
       exit(1)
   
 logging.info('Export at %s ', foreground_folder )
else:
 logging.info('Skip steps 4 and 5: filtering_result_by_expression, plot_heatmaps ')

if runRemoveTempFile:
 #6. clean temporary files
 logging.info("")
 logging.info("Step 6: Remove temporary files in background or foreground calculations. ")

 result_code=os.system('bpb3 clean_tmp --background_out_file ' + background_folder + 
                ' --foreground_out_file ' + foreground_folder  +  ' --chemical_potentials ' + chemical_potentials  +
                ' --background_pwm_file ' + background_pwm_folder + ' --clean_both 2' )
 if result_code !=0:
   logging.info("Error in clean tmp ")
   exit(1)
else:
 logging.info('Skip step 6: remove temporary file %s ', background_pwm_folder) 













