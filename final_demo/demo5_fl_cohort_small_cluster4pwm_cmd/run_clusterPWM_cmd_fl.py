import os
import glob
import sys
from bpb3.script.script_high.other import common
import time
#this script is used to evaluate effect of TF binding affinity changes due to DNA sequence mutations based on clustered PWMs from abc4pwm - first level analysis
#Here, most of bpb3 command line functions and parameter configures are shown in the same file
#
########Project main path parameters ##############
#Project Output path
out_data_folder = '../../final_demo_data/demo5_fl_cohort_small_cluster4pwm_cmd/out/'
common.check_folder(out_data_folder)
patient_data_folder = os.path.join(out_data_folder, '../',"patient_data")
normal_rnaseq_folder = os.path.join(out_data_folder,'../', "normal_gcb_counts") 
project_name= 'demo5_fl_cohort_small_cluster4pwm_cmd'

#find patient muation files
patient_mut_files = []
for pf in glob.glob(os.path.join(patient_data_folder, "D*/")):
      patient_mut_files.append(os.path.join(pf, "icgc_mut*.tsv"))

#Genome path
genome_folder=os.path.abspath("../../final_demo_data/genome_data/hg19")
genome_fasta_file=  os.path.join(genome_folder, "human_g1k_v37_decoy.fasta") 
gene_annotation_file= os.path.join(genome_folder, "gencode.v19.annotation.gtf")
gene_length_file =  os.path.join(genome_folder, "gene_lengths.tsv")


# output folders for mussd
mussd_result_folder=os.path.join(out_data_folder, 'mussd_blocks')
foreground_folder = os.path.join(out_data_folder, "foreground")
background_folder = os.path.join(out_data_folder, "background")
background_pwm_folder = os.path.join(out_data_folder, "pwm_for_background")
 
#Input gene expression folder and files
normal_count_files =  glob.glob(os.path.join(normal_rnaseq_folder, "*.only_counts.tsv"))
if len(normal_count_files)<1:
    print('Normal count file, ', normal_rnaseq_folder, ' not find!')
    exit(1)
donor_count_files = glob.glob(os.path.join(patient_data_folder, "D*/D*expression.tsv"))
if len(donor_count_files)<1:
    print('Tumor count file, ', donor_count_files, ' not find!') 
    exit(1)

#Setup logger file
logging =common.setup_logger(out_data_folder,project_name, time.time())
logging.info("")

# make cluster pwms 
#  Generate tempary pwms based on clustered and uncertain clustered PWMs ############
in_clustered_pwms_path=os.path.join(out_data_folder,'../')
in_clustered_pwms_folder_name='abc4pwm_quality_assessed_out_seed5'
in_uncertain_pwms_folder_name='abc4pwm_uncertain_pwms'
out_tmp_pwms_folder_name='tmp_pwm'
top_ranked_TF_from_clusteredPWM= 15

# folder with PWM models for TF binding affinity
#pwm_folder = os.path.join("../../public_data/pwm")
pwm_folder=os.path.join(in_clustered_pwms_path,out_tmp_pwms_folder_name)

#0. prepare clustered PWMs file for bpb3
cmd='bpb3 make_cluster4pwm --in_pwm_path ' + in_clustered_pwms_path +  \
        ' --in_clustered_pwm_folder ' +  in_clustered_pwms_folder_name + \
        ' --out_file_folder ' + out_tmp_pwms_folder_name + \
        ' --in_uncertain_pwm_folder ' + in_uncertain_pwms_folder_name

result_code=os.system(cmd)
if result_code != 0:
      print("bpb3 make_cluster4pwm failed. See message in the console to find what went wrong.")
      exit(1)
logging.info('Export at %s', os.path.join(in_clustered_pwms_path, out_tmp_pwms_folder_name))
pwm_folder=os.path.join(out_data_folder,'../',out_tmp_pwms_folder_name)


##############BPB3 pipeline parameters in each STEP ############################
#1. Differential gene expression options
differentially_expressed_genes_file = os.path.join(out_data_folder, "differentially_expressed_genes.txt")
group1_str='tumor'
group2_str='normal'
#Export file path
group1_median_rpkm_file = os.path.join(out_data_folder, group1_str + "_median_patient_rpkm.tsv")
group2_median_rpkm_file = os.path.join(out_data_folder, group2_str +"_median_patient_rpkm.tsv")

# maximum Kolmogorov-Smirnov test P value to consider genes as differentially expressed
differential_expression_p = 0.05
# whether to do quantile normalization of RPKM values, True or False
differential_expression_quantile_normalization = True
# whether to do log transformation of RPKM values, True or False
differential_expression_log_transform = True
# whether to do z-score transformation of RPKM values, True or False
differential_expression_z_score_transform = False
# minimum fold change of [quantile normalized] RPKM to consider genes as differentially expressed. None to disable
differential_expression_min_fold_change = None
#export full expression data
isExportAll = True
#test_method 1 for T-test
test_method=1
#minimum median RPKM in each group
minimum_median_RPKM=0  #1.0

# 2. Region of interest selection options
# number of base pairs to take upstream of TSS
upstream_size = 1000
# number of base pairs to take downstream of TSS
downstream_size = 1000
#Export file path or selected enhancer regions file
regions_file = os.path.join(out_data_folder, "regions.bed")

# Here, we do not use bpb3 function to extract a specified TSS regions, but manually input a bed format region file
# for DNA muatation enrichment and further TF binding affinity change significance test 
#regions_file='../../public_data/mr_test/in/mr_intersect_enhancer_regions_chr18noChr.tsv'

# 3. MuSSD options
#in this example, we include all patiaints and all available mutations in the regions
# minimum number of patients in a hot mutation region
min_patients_in_block = 3
# minimum number of mutations in a hot mutation region
min_block_size = 3
# maximum distance between mutations in a hot mutation region
cluster_distance = 30
block_distance=500
# number of base pairs to add to the left and right side of a mutation block when extracting sequences
block_flank = 25

#4. P values for significant muation blocks
p_value_for_significant_blocks=0.001
mussd_command_file=os.path.join(out_data_folder, 'run_mussd.sh')
pval_correction=True
significant_blocks_file = os.path.join(out_data_folder, "significant_blocks.txt")

# 5. BayesPI-BAR options
# chemical potentials to use
chemical_potentials = "none -10 -13 -15 -18 -20"
# sequence shuffling iterations for dbA calculation
shuffling_iterations = 10000
#Parallel option path
parallel_options_file =  os.path.abspath(os.path.join(os.path.dirname(__file__), "parallel_options.txt"))

# 6. Background model options
# max rank of a PWM file in the foreground results to include in the background model
max_ranking_for_background = 30
resampling_iterations = 10
gene_samples = 8
background_shuffling_iterations = 1000
background_command_file = os.path.join(out_data_folder, "make_background.sh")

# 7. mutation signature definition for generating background mutations, set to None to disable
mutation_signature_file = None  # "../../data/signature_7.tsv"
# set to True to use tumor mutations from patients as background mutations, False otherwise.
# This and the mutation signature are mutually exclusive.
# In both are disabled, uniform mutations will be generated
tumor_mutations_as_background = True

# Results selection options for significant affinity changes
# P value for ranksum test in the output
affinity_change_p_value = 0.05
pval_correction4affinity=True

#8. True or False parameters to run the bpb3, by either enable  or disable one of bpb3 steps
runDiffExp= True
runRegionsFile=True
runMussd= True
runHighMutBlocks= True
runBayesPiBAR= True
runBackgroundModel= True
runBlockSignificance= True

filterTFbyGene=True
runRemoveTempFile= True


############Usually, below lines do not need to be changed in order to run it on a new data set ############## 

################START CALC ################
# 1. Differential gene expression
if runDiffExp:
 logging.info("Step 1: Computing differentially expressed genes")
 result_code = os.system("bpb3 differential_expression " +
                                " --group_2_count_files " + " ".join(normal_count_files) +
                                " --group_1_count_files " + " ".join(donor_count_files) +
                                " --gene_lengths " + gene_length_file +
                                " --output_file " + differentially_expressed_genes_file +
                                " --output_group_1_rpkm " + group1_median_rpkm_file +
                                " --output_group_2_rpkm " + group2_median_rpkm_file +
                                " --p_value " + str(differential_expression_p) +
                                " --test_method " + str(test_method) +
                                " --quantile_normalization" * differential_expression_quantile_normalization +
                                " --log_transform" * differential_expression_log_transform +
                                " --z_score" * differential_expression_z_score_transform +
                                " --output_all_values " * isExportAll + 
                                (" --min_medianRPKM " + str(minimum_median_RPKM) 
                                 if minimum_median_RPKM is not None else "") + 
                                (" --min_fold_change " + str(differential_expression_min_fold_change)
                                 if differential_expression_min_fold_change is not None else ""))
 if result_code !=0:
       logging.info("Error in differential expresssion ")
       exit(1)
  
 logging.info('Export at %s ', differentially_expressed_genes_file )
else :
 logging.info('Skip step 1: differential gene expression analysis %s ', differentially_expressed_genes_file)
 

#2. Make regions of interest near promoters of differentially expressed genes
if runRegionsFile: 
  cmd0="bpb3 gene_regions " +\
          " --gene_annotation " + gene_annotation_file +\
          " --transform_chrom" \
          " --output_file " + regions_file + \
          " --selected_genes " + differentially_expressed_genes_file + \
          " --upstream_size " + str(upstream_size) + \
          " --downstream_size " + str(downstream_size)
  result_code=os.system(cmd0)
  if result_code !=0:
       logging.info("Error in bpb3 gene_regions ")
       exit(1)
  else:
        logging.info('Export at %s ', regions_file )
else:
  logging.info('Skip step 2 regions export, but use selected regions for MUSSD %s  ', regions_file)

if runMussd:
 # 3. Run MuSSD
 logging.info("")
 logging.info("Step 3: Calculating mutation blocks with MuSSD")
 if not os.path.exists(mussd_result_folder):
   cmd_head='#!/bin/bash' 
   cmd ="bpb3 mussd --patient_mutations " + " ".join(patient_mut_files) + " \\\n\t" + \
                   " --result_folder " + mussd_result_folder +  " \\\n\t" + \
                   " --genome " + genome_fasta_file +  " \\\n\t" + \
                   " --min_patients_in_block " + str(min_patients_in_block) +  " \\\n\t" + \
                   " --block_flank " + str(block_flank) +  " \\\n\t" + \
                   " --min_block_size " + str(min_block_size) +  " \\\n\t" + \
                   " --cluster_distance " + str(cluster_distance) +  " \\\n\t" + \
                   " --block_distance " + str(block_distance) + "\\\n\t" + \
                   " --regions " + regions_file
   with open(mussd_command_file,'w') as f:
       f.write(cmd_head + "\n" + cmd + "\n")
   f.close()
   #if input VCF files too many then have to run it in bash 
   result_code= os.system("bash "+ mussd_command_file)
   if result_code !=0:
       #print(cmd)
       logging.info("Error in mussd ")
       exit(1)
  
   logging.info('Export at %s ', mussd_result_folder)
   logging.info("Information about patients is in the file %s", os.path.join(mussd_result_folder,
                                                                                  "patients_summary.tsv"))
   logging.info("Information about mutations is in the file %s", os.path.join(mussd_result_folder,
                                                                                   "mutations_summary.tsv"))
   logging.info("Information about found mutation blocks is in the file %s", os.path.join(mussd_result_folder,
                                                                                            "blocks_summary.tsv"))
 else:
   logging.info("MuSSD result is already present, skipping step 3")

 logging.info("There are %d patients in the mutation dataset", len(patient_mut_files))
else:
 logging.info('Skip step 3 MuSSD %s ', mussd_result_folder)

if runHighMutBlocks:
 # 4. Select blocks that are mutated significantly more than average
 #significant_blocks_file = os.path.join(out_data_folder, "significant_blocks4mr_enhancer_regions.txt")
 logging.info("")
 logging.info("Step 4: Selecting blocks that are mutated more than average")
 result_code = os.system("bpb3 highly_mutated_blocks " +
                            " --blocks_folder " + mussd_result_folder +
                            " --p_value " + str(p_value_for_significant_blocks) + 
                            " --output_file " + significant_blocks_file +
                            " --correctPval " * pval_correction ) 
 if result_code !=0:
       logging.info("Error in highly mutated blocks ")
       exit(1)

 logging.info('Export at %s ', significant_blocks_file )
else:
 logging.info('Skip step 4 highly_mutated_blocks test %s ', significant_blocks_file )


if os.stat(significant_blocks_file).st_size ==0:
  logging.info('No significant mutation blocks find, in P< %s ', str(p_value_for_siginificant_blocks) )
  logging.info('Program stop!!')
  exit(1)


#Read  significant mutation blocks
significant_block_tags = frozenset(common.read_lines(significant_blocks_file))
logging.info("%d blocks have significantly more mutations than expected: %s", len(significant_block_tags),
                 ", ".join(significant_block_tags))

#if runBayesPiBAR:
 # 5. Compute foreground for each significant block
logging.info("")
logging.info("Step 5: Running BayesPI-BAR on patient-specific mutation blocks")

common.check_folder(foreground_folder)
block_res_folders=[]
for i, block_tag in enumerate(significant_block_tags):
        block_res_folder = os.path.join(foreground_folder, block_tag)
        ref_seq = os.path.join(mussd_result_folder, block_tag + "_ref.fasta")
        alt_seq = os.path.join(mussd_result_folder, block_tag + "_alt.fasta")
        block_res_folders.append(block_res_folder)
        if os.path.exists(block_res_folder):
            logging.info("The foreground output for block %s [%d of %d] already exists. BayesPI-BAR will now check its "
                         "consistency and recompute the broken or missing data ", block_tag, i + 1,
                         len(significant_block_tags))
        else:
            logging.info("BayesPI-BAR will now compute the delta-dbA scores and rankings for block %s [%d of %d]",
                         block_tag, i + 1, len(significant_block_tags))
        if runBayesPiBAR:
           result_code = os.system("bpb3 bayespi_bar " +
                                " --reference_sequences " + ref_seq +
                                " --alternate_sequences " + alt_seq +
                                " --pwm_folder " + pwm_folder +
                                " --chemical_potentials " + chemical_potentials +
                                " --iterations " + str(shuffling_iterations) +
                                " --result_folder " + block_res_folder +
                                " --max_rank " + str(max_ranking_for_background) +
                                " --reuse_output" +
                                " @" + parallel_options_file)

           if result_code !=0:
              logging.info("Error in bayespi bar ")
              exit(1)
        
           logging.info('Export at %s', block_res_folder )
        else:
           logging.info('Skip step 5 bayePI-BAR %s ' , block_res_folder) 

if runBackgroundModel:
 # 6. Compute the background
 logging.info("")
 logging.info("Step 6: Computing the background model")
 logging.info("Creating shell file with the command to compute the background: %s", background_command_file)
 result_code = os.system("bpb3 choose_background_parameters " +
                            " --blocks_folder " + mussd_result_folder +
                            " --foreground_folder " + foreground_folder +
                            " --background_pwm_folder " + background_pwm_folder +
                            " --output_file " + background_command_file +
                            " --block_ids " + " ".join(significant_block_tags) +
                            " --genome " + genome_fasta_file +
                            " --regions " + regions_file +
                            " --pwm_folder " + pwm_folder +
                            " --iterations " + str(background_shuffling_iterations) +
                            " --chemical_potentials " + chemical_potentials +
                            " --background_output_folder " + background_folder +
                            " --block_resample_iterations " + str(resampling_iterations) +
                            " --block_samples_to_take " + str(gene_samples) +
                            " --reuse_output" +
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
 logging.info('Skip step 6 background_computation %s ', background_folder)


if runBlockSignificance:
 # 7. Run significance tests
 logging.info("")
 logging.info("Step 7: Significance tests will now be computed comparing background and foreground delta-dbA "
                 "distributions for each PWM in each mutation block")

 for block_tag in significant_block_tags:
        block_folder = os.path.join(foreground_folder, block_tag)
        result_file = os.path.join(block_folder, "result.tsv")
        print( block_tag, ":" )
        result_code = os.system("bpb3 affinity_change_significance_test " +
                                " --background_folder " + os.path.join(background_folder, "bayespi_bar_result") +
                                " --foreground_folder " + block_folder +
                                " --output_file " + result_file +
                                " --max_rank " + str(max_ranking_for_background) +
                                " --p_value " + str(affinity_change_p_value) +
                                " --pval_correction " * pval_correction4affinity)
        if result_code !=0:
            logging.info("Error in affinity change significance test %s", block_tag)
            exit(1)
        
        logging.info("The table of significantly affected PWMs for block %s is stored in the file %s ", block_tag,
                     result_file)

if filterTFbyGene:
 # 8. Filter results by TF gene expression
 logging.info("")
 logging.info("Step 8: PWMs corresponding to the non-expressed TFs will now be removed from the results")

 for block_tag in significant_block_tags:
        block_folder = os.path.join(foreground_folder, block_tag)
        result_file = os.path.join(block_folder, "result.tsv")
        filtered_result_file = os.path.join(block_folder, "result_filtered.tsv")

        print("Filtering results by gene expression for block", block_tag, "...")
        result_code = os.system("bpb3 filter_results_by_gene_expression " +
                                " --results_file " + result_file +
                                " --expression " + group1_median_rpkm_file +
                                " --output_file " + filtered_result_file)

        if result_code != 0:
            logging.info("Filtering results for block %s failed. "
                         "See messages in the console to find what went wrong.", block_tag)
            exit(1)

        logging.info("The table of significantly affected PWMs for block %s, containing only those TFs which are "
                     "expressed in patients, is stored in the file %s", block_tag,
                     filtered_result_file)

 #9. plot all heatmaps
 logging.info("")
 logging.info("Step 9: Generate plots or heatmaps for the results")

 result_code=os.system('bpb3 make_plots --foreground_folder ' + foreground_folder +
          ' --background_folder ' + background_folder + 
          ' --mussd_result_folder ' + mussd_result_folder )
 #bpb3 make_plots --foreground_folder '../../data/skin_cancer/out/foreground'  --background_folder '../../data/skin_cancer/out/pwm_for_background' --mussd_result_folder '../../data/skin_cancer/out/mussd_blocks'

 if result_code !=0:
       logging.info("Error in make plots ")
       exit(1)
   
 logging.info('Export at %s ', foreground_folder )
else:
 logging.info('Skip steps 7, 8, and 9 affinity_significance_test, filtering_result_by_expression, plot_heatmaps ')

if runRemoveTempFile:
 #10. clean temporary files
 logging.info("")
 logging.info("Step 10: Remove temporary files in both background and foreground calculations. ")

 result_code=os.system('bpb3 clean_tmp --background_out_file ' + background_folder + 
                ' --foreground_out_file ' + foreground_folder  +  ' --chemical_potentials ' + chemical_potentials + 
                ' --background_pwm_file ' + background_pwm_folder + ' --clean_both 2' )
 if result_code !=0:
   logging.info("Error in clean tmp ")
   exit(1)
else:
 logging.info('Skip step 10 remove temporary file %s ', background_pwm_folder) 


