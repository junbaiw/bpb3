import os
import glob
from bpb3.script.script_high.other import common
import time

#this script is used to evaluate mutations in DNA methylation regions (MRs) that are overlapping with enhancers!!
#Here, most of bpb3 command line functions and parameter configures are shown in the same file

###########Setup main Path ################
#Project Output path
out_data_folder = os.path.abspath("../../final_demo_data/demo7_mr_enhancers/out")

#input data path
patient_data_folder = os.path.join(out_data_folder,'../', "patient_data")
normal_rnaseq_folder = os.path.join(out_data_folder,'../' ,"normal_gcb_counts")

#project name
project_name='demo7_mr_enhancers'

#Genome path
#genome data too large which needs to be stored in a separeted place from the demo data in github
genome_folder=os.path.abspath("../../final_demo_data/genome_data/hg19/")
genome_fasta_file= os.path.join(genome_folder, "human_g1k_v37_decoy.fasta")
gene_annotation_file= os.path.join(genome_folder, "gencode.v19.annotation.gtf")

#Start to make folder for clustered pwms, which is the only difference if input clustered PWMs for bayesPI-BAR ! 
#this part works when runCluster4PWM is True 
#  Generate tempary pwms based on clustered and uncertain clustered PWMs ############
in_clustered_pwms_path=os.path.join(out_data_folder,'../')
in_clustered_pwms_folder_name='abc4pwm_quality_assessed_out_seed5'
in_uncertain_pwms_folder_name='abc4pwm_uncertain_pwms'
out_tmp_pwms_folder_name='tmp_pwm'
#end custer4PWM parameters 

# folder with clustered PWM models for TF binding affinity
#pwm_folder=os.path.join(in_clustered_pwms_path,out_tmp_pwms_folder_name)
# folder with original PWM models for TF binding affinity
pwm_folder = os.path.join("../../final_demo_data/pwm")

# output data folders
foreground_folder = os.path.join(out_data_folder, "foreground")
background_folder = os.path.join(out_data_folder, "background")
background_pwm_folder = os.path.join(out_data_folder, "pwm_for_background")
mussd_result_folder = os.path.join(out_data_folder, "mussd_blocks")

#normal gene expression file name under folder normal_rnaseq_folder
normal_count_files_postprefix="*.only_counts.tsv"
gene_length_file = os.path.join(genome_folder, "gene_lengths.tsv")
#gene expression file name postprefix in folder patient_data_folder
donor_count_files_postprefix="D*/*_expression.tsv"

###########Setup Pipeline parameters #########
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
minimum_median_RPKM= 0   #1.0

# 2. Region of interest selection options
# number of base pairs to take upstream of TSS
#upstream_size = 1000
# number of base pairs to take downstream of TSS
#downstream_size = 1000
#Export file path
#regions_file = os.path.join(out_data_folder, "regions.bed")
#region selection function is disabled , a predefined region is used such as DMR overlapping to enhancers with a 5kb flank region in both side 
regions_file='../../final_demo_data/demo7_mr_enhancers/in/dmr_intersect_enhancer_regions_chr18noChr_5kb_flank.tsv'


# 3. MuSSD options
# minimum number of patients in a hot mutation region
min_patients_in_block = 3
# minimum number of mutations in a hot mutation region
min_block_size = 3
# maximum distance between mutations in a hot mutation region
cluster_distance = 30
block_distance=500
# number of base pairs to add to the left and right side of a mutation block when extracting sequences
block_flank = 25

#4. P values for significant muation blocks use a big positive value >1 to switch off this filtering
p_value_for_significant_blocks=0.001
mussd_command_file=os.path.join(out_data_folder, 'run_mussd.sh')
pval_correction= True

#patient mutation file names under folder patient_data_folder
patient_mut_files_postprefix="DO*/icgc_mutations*.tsv"

#Export file path of significant mutation blocks
significant_blocks_file = os.path.join(out_data_folder, "significant_blocks.txt")

# 5. BayesPI-BAR options
# chemical potentials to use
chemical_potentials = "none -10 -13 -15 -18 -20"
# sequence shuffling iterations for dbA calculation
shuffling_iterations = 10000
#Parallel option path
parallel_options_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "parallel_options.txt"))

# 6. Background model options
# max rank of a PWM file in the foreground results to include in the background model
max_ranking_for_background = 10
resampling_iterations = 10
gene_samples = 100
background_shuffling_iterations = 1000
background_command_file = os.path.join(out_data_folder, "make_background.sh")

# 7. mutation signature definition for generating background mutations, set to None to disable
mutation_signature_file = None  # "../../final_demo_data/background_data/signature_7.tsv"
# set to True to use tumor mutations from patients as background mutations, False otherwise.
# This and the mutation signature are mutually exclusive.
# In both are disabled, uniform mutations will be generated
tumor_mutations_as_background = True

# Results selection options for significant affinity changes
# P value for ranksum test in the output
affinity_change_p_value = 0.05
pval_correction4affinity= True
max_ranking_for_significant_blocks=30

#8. True or False parameters to run the bpb3, by either enable  or disable one of bpb3 steps

#step 1. differential expression analysis
runDiffExp= True

#step 2: region extraction, disable this function becaues we use predefined DMR overlapping enhancer regions
runRegionsFile=False


#step 3: mussd mutation block detection
runMussd= True

#step 4: highly mutated mutation block test
runHighMutBlocks= True

#step 5: do foreground bayespi-bar calculation
runBayesPiBAR= True

#step 6: do background bayespi-bar calcuation
runBackgroundModel= True

#step 7: do TF affinity change significance test
runBlockSignificance= True

#step 8: filter TF by gene expression level
filterTFbyGene= True

#step 9: plot heatmap for TFs with significaly changed binding affinity
isPlot=True

#step 10: remove temporary files in both foreground and background calculations
runRemoveTempFile=True

#option step: if the input PWMs are clustered PWMs from abc4pwm, then this option shall be True
runCluster4PWM=False



