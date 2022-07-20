#This bpb3 parameter configure is for demo purpose
import os
import glob
from bpb3.script.script_high.other import common
import time

###########Setup main Paths and folders ################
#Project main Output path
out_data_folder = os.path.abspath("../../final_demo_data/demo4_fl_cohort_small_cluster4pwm/out")

#Project input path for patient mutation data and gene expression data
patient_data_folder = os.path.join(out_data_folder,'../', "patient_data")

#Project input path for gene expression data of normal control samples
normal_rnaseq_folder = os.path.join(out_data_folder,'../', "normal_gcb_counts")

#Project name in logo file
project_name='demo4_fl_cohort_smallcluster4pwm'

#Input Genome path
genome_folder=os.path.abspath("../../final_demo_data/genome_data/hg19/")

#Input human genome sequence
genome_fasta_file= os.path.join(genome_folder, "human_g1k_v37_decoy.fasta")

#Input gtf file of human genes from GenCode
gene_annotation_file= os.path.join(genome_folder, "gencode.v19.annotation.gtf")

#Project output folders for foreground and background calculations
foreground_folder = os.path.join(out_data_folder, "foreground")
background_folder = os.path.join(out_data_folder, "background")

#Background pwm folder 
background_pwm_folder = os.path.join(out_data_folder, "pwm_for_background")

#output folder of predicted muation blocks by MUSSD algorithm
mussd_result_folder = os.path.join(out_data_folder, "mussd_blocks")

##### Set up input file name postprefix #########
#Input gene expression folder is located in patient data folder
#postprefix of gene expression file name for normal control samples such as files in folder normal_rnaseq_folder
normal_count_files_postprefix="*.only_counts.tsv"

#gene length file name
gene_length_file = os.path.join(genome_folder, "gene_lengths.tsv")

#postprefix of gene expression file name for patients such as files in folder patient_data_folder
donor_count_files_postprefix="D*/*_expression.tsv"

#### 0. set up PWM file path for clustered PWMs or common PWMs 
##################Here is the only difference if input PWMs are either clustered PWMs or common PWMs! ################ 
################# Generate tempary pwms based on clustered and uncertain clustered PWMs ############
#Input path of clustered PWMs
in_clustered_pwms_path=os.path.join(out_data_folder,'../')

# file folder of clustered PWMs passed quality control in abc4pwm package 
in_clustered_pwms_folder_name='abc4pwm_quality_assessed_out_seed5'

# file folder of PWMs do not pass quality control in abc4pwm
in_uncertain_pwms_folder_name='abc4pwm_uncertain_pwms'

# new input folder of clustered PWMs when reorganizing PWMs file names based on the two aforementioned folders.
# these renamed PWMs will be used as input PWMs for applying clustered  PWMs in bpb3 
out_tmp_pwms_folder_name='tmp_pwm'

# final pwm folder path for TF binding affinity calculation.
# if we do not use clustered PWMs as input, then this can be simply replaced by a folder path of common PWMs such as the folder contains 1772 PWMs in bpb3
pwm_folder=os.path.join(in_clustered_pwms_path,out_tmp_pwms_folder_name)

###########Setup Pipeline parameters from step 1 to step 8 before running bpb3 #########
#### 1. Differential gene expression options
#export file name for differentially expressed genes
differentially_expressed_genes_file = os.path.join(out_data_folder, "differentially_expressed_genes.txt")

#group name for tumor and normal samples
group1_str='tumor'
group2_str='normal'
#Export gene expression files path
group1_median_rpkm_file = os.path.join(out_data_folder, group1_str + "_median_patient_rpkm.tsv")
group2_median_rpkm_file = os.path.join(out_data_folder, group2_str +"_median_patient_rpkm.tsv")

# maximum Kolmogorov-Smirnov test or other statistic test P value to consider genes as differentially expressed
differential_expression_p = 0.05
# whether to do quantile normalization of RPKM values, True or False
differential_expression_quantile_normalization = True
# whether to do log transformation of RPKM values, True or False
differential_expression_log_transform = True
# whether to do z-score transformation of RPKM values, True or False
differential_expression_z_score_transform = False
# minimum fold change of [quantile normalized] RPKM to consider genes as differentially expressed. None to disable
differential_expression_min_fold_change = None
#export full expression data , True or False
isExportAll = True
#select test method:  1 for T-test, 0 for KS-test
test_method=1
#minimum median RPKM allowed in each group, input 0 to disable this option
minimum_median_RPKM=0  #1.0

#### 2. Region of interest selection options
# number of base pairs to take upstream of TSS
upstream_size = 1000
# number of base pairs to take downstream of TSS
downstream_size = 1000
#Export file path for extracted regions
regions_file = os.path.join(out_data_folder, "regions.bed")

# If we do not use bpb3 function to extract a specified TSS regions, then manually input a bed format region file
# for example, using enhancer regions overlapping with Differential Methylation Regions for DNA muatation enrichment and the TF binding affinity change significance test 
#regions_file='../demo7_mr_enhancers/in/dmr_intersect_enhancer_regions_chr18noChr_5kb_flank.tsv'

#### 3. MuSSD options
# minimum number of patients in a hot mutation region
min_patients_in_block = 3
# minimum number of mutations in a hot mutation region
min_block_size = 3 
# maximum distance between mutations in a hot mutation region
cluster_distance = 30
# maximum distance between two adjacent mutaiton blocks
block_distance=500
# number of base pairs to add to the left and right side of a mutation block when extracting sequences
block_flank = 25

#### 4. P values for significant muation blocks
p_value_for_significant_blocks=0.001
#commamnd line file name for running mussd
mussd_command_file=os.path.join(out_data_folder, 'run_mussd.sh')
#whether to use bonferroni corrected p value, True or False
pval_correction=True

#patient mutation file names under folder patient_data_folder
patient_mut_files_postprefix="DO*/icgc_mutations*.tsv"

#Export file path of significant mutation blocks
significant_blocks_file = os.path.join(out_data_folder, "significant_blocks.txt")

#### 5. BayesPI-BAR options
# a range of chemical potentials to use
chemical_potentials = "none -10 -13 -15 -18 -20"
# sequence shuffling iterations for dbA calculation
shuffling_iterations = 10000
#Parallel option path
parallel_options_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "parallel_options.txt"))
#maximum number of TFs will be ranked in foreground calculation
max_ranking_for_foreground=30
#whether both forground and background ddba are already exists, True or False
start_from_integration=False

#### 6. Background model options
# max rank of a PWM file in the foreground results to include in the background model
max_ranking_for_background = 20
# number of random resampling
resampling_iterations = 10
# number of genes will be sampled
gene_samples = 8
# number of random shuffling in backgroun calculation
background_shuffling_iterations = 1000
# command line file for running background calculation
background_command_file = os.path.join(out_data_folder, "make_background.sh")

#### 7. mutation signature definition for generating background mutations, set to None to disable
mutation_signature_file = None  # "../../data/signature_7.tsv"
# set to True to use tumor mutations from patients as background mutations, False otherwise.
# This and the mutation signature are mutually exclusive.
# In both are disabled, uniform mutations will be generated
tumor_mutations_as_background = True

# Results selection options for significant affinity changes
# P value for ranksum test in the output
affinity_change_p_value = 0.05
# whether to use bonferroni corrected p value, True or False
pval_correction4affinity=True
#signifcant block parameters, maximum number of top ranked PWMs
max_ranking_for_significant_blocks= 30

###### End of Parameter setup, below options for running each of the analysis pipeline in bpb3 ############
#### 8. Use True or False to either enable or disable the corresponding analysis step when running bpb3 #####

#step 1. differential expression analysis
runDiffExp= True

#step 2: region extraction
runRegionsFile= True

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
filterTFbyGene=True

#step 9: plot heatmap for TFs with significaly changed binding affinity
isPlot=True

#step 10: remove temporary files in both foreground and background calculations
runRemoveTempFile= True

#option step: if the input PWMs are clustered PWMs from abc4pwm, then this option shall be True, otherwise it is False.
runCluster4PWM=True


