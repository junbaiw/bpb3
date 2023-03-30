import os
import glob
import sys
from bpb3.script.script_high.other import common
import time
import argparse

def my_parser(parser):
  required= parser.add_argument_group("Required ")
  required.add_argument('--import_bpb3_config', help="File for import bpb3 configure file , and the name shall be bpb3_config.py ",
                        required= True, type=argparse.FileType(), metavar="FILE")
  return parser

def run(args):
 
 #read configure file name and its path
 configure_path=os.path.split(os.path.abspath(args.import_bpb3_config.name))
 configure_file= configure_path[1].replace('.py','')
 print('Add path \n', configure_path[0])
 sys.path.insert(0, configure_path[0])

 #load bpb3 configure file
 #this works but input configure file name cant be changed. 
 #import bpb3_config

 #this can use any name for input configure file
 cmd='''\"import ''' + configure_file + ''' as bpb3_config\"'''
 print(cmd)
 eval("exec(" + cmd +", globals(),globals() )")

 ###Load Parameters##########
 #Usually, this script does not need change because all parameers are setup in bpb3_config.py

 ###########Setup main Path ################
 #Project Output path
 out_data_folder = bpb3_config.out_data_folder 
 patient_data_folder = bpb3_config.patient_data_folder
 normal_rnaseq_folder = bpb3_config.normal_rnaseq_folder
 project_name= bpb3_config.project_name

 #Setup logger file
 logging =common.setup_logger(out_data_folder,project_name, time.time())

 #Genome path
 genome_folder=bpb3_config.genome_folder
 genome_fasta_file= bpb3_config.genome_fasta_file
 gene_annotation_file= bpb3_config.gene_annotation_file

 # folder with PWM models for TF binding affinity
 #pwm_folder = bpb3_config.pwm_folder

 # output folders
 foreground_folder = bpb3_config.foreground_folder
 background_folder = bpb3_config.background_folder
 background_pwm_folder = bpb3_config.background_pwm_folder
 
 mussd_result_folder = bpb3_config.mussd_result_folder

 #Input gene expression folder
 normal_count_files =glob.glob(os.path.join(normal_rnaseq_folder, bpb3_config.normal_count_files_postprefix))
 #normal_count_files= bpb3_config.get_normal_count_files()
 gene_length_file = bpb3_config.gene_length_file
 #here may need change!!
 #donor_count_files = bpb3_config.get_donor_count_files()
 donor_count_files=glob.glob(os.path.join(patient_data_folder,bpb3_config.donor_count_files_postprefix))

 ###########Setup Pipeline parameters #########33
 #1. Differential gene expression options
 differentially_expressed_genes_file = bpb3_config.differentially_expressed_genes_file 
 #group1_str='tumor'
 #group2_str='normal'
 #Export file path
 group1_median_rpkm_file = bpb3_config.group1_median_rpkm_file
 group2_median_rpkm_file =  bpb3_config.group2_median_rpkm_file

 # maximum Kolmogorov-Smirnov test P value to consider genes as differentially expressed
 differential_expression_p = bpb3_config.differential_expression_p

 # whether to do quantile normalization of RPKM values, True or False
 differential_expression_quantile_normalization =  bpb3_config.differential_expression_quantile_normalization

 # whether to do log transformation of RPKM values, True or False
 differential_expression_log_transform =  bpb3_config.differential_expression_log_transform

 # whether to do z-score transformation of RPKM values, True or False
 differential_expression_z_score_transform =  bpb3_config.differential_expression_z_score_transform

 # minimum fold change of [quantile normalized] RPKM to consider genes as differentially expressed. None to disable
 differential_expression_min_fold_change = bpb3_config.differential_expression_min_fold_change
 #export full expression data
 isExportAll =  bpb3_config.isExportAll
 #test_method 1 for T-test
 test_method=  bpb3_config.test_method
 #minimum median RPKM in each group
 minimum_median_RPKM= bpb3_config.minimum_median_RPKM

 # 2. Region of interest selection options
 # number of base pairs to take upstream of TSS
 if hasattr(bpb3_config,'upstream_size'):
   upstream_size =  bpb3_config.upstream_size
 else:
   upstream_size=1000

 # number of base pairs to take downstream of TSS
 if hasattr(bpb3_config,'downstream_size'):
   downstream_size =  bpb3_config.downstream_size
 else:
   downstream_size= 1000
 #Export file path
 regions_file =  bpb3_config.regions_file

 # 3. MuSSD options
 # minimum number of patients in a hot mutation region
 min_patients_in_block = bpb3_config.min_patients_in_block
 # minimum number of mutations in a hot mutation region
 min_block_size = bpb3_config.min_block_size 
 # maximum distance between mutations in a hot mutation region
 cluster_distance = bpb3_config.cluster_distance
 # number of base pairs to add to the left and right side of a mutation block when extracting sequences
 block_flank = bpb3_config.block_flank
 if hasattr(bpb3_config, 'block_distance'):
    block_distance=bpb3_config.block_distance
 else:
    #1000 is default setting in mussd
    block_distance=1000

 #4. P values for significant muation blocks
 p_value_for_significant_blocks=  bpb3_config.p_value_for_significant_blocks
 if hasattr(bpb3_config,'pval_correction'):
   pval_correction= bpb3_config.pval_correction
   if pval_correction:
      logging.info('Assign P-value correction in significant mutation blocks test')
   else:
      logging.info('No P-value correction in significant mutation blocks test')
 else:
   logging.info('Default False p-value correction in significant mutation blocks test')
   pval_correction=False
 mussd_command_file=  bpb3_config.mussd_command_file

 #patient_mut_files = bpb3_config.patient_mut_files
 patient_mut_files=glob.glob(os.path.join(patient_data_folder,bpb3_config.patient_mut_files_postprefix))
  

 #Export file path of significant mutation blocks
 significant_blocks_file = bpb3_config.significant_blocks_file  

 # 5. BayesPI-BAR options
 # chemical potentials to use
 chemical_potentials = bpb3_config.chemical_potentials
 # sequence shuffling iterations for dbA calculation
 shuffling_iterations = bpb3_config.shuffling_iterations
 #Parallel option path
 parallel_options_file = bpb3_config.parallel_options_file
 #ranking in bayespi-bar 
 if hasattr(bpb3_config, 'max_ranking_for_foreground'):
   max_ranking_for_foreground= bpb3_config.max_ranking_for_foreground
 else:
   max_ranking_for_foreground=10
 #is skip ddba calculation
 if hasattr(bpb3_config, 'start_from_integration'):
   start_from_integration=bpb3_config.start_from_integration
 else:
   start_from_integration=False

 # 6. Background model options
 # max rank of a PWM file in the foreground results to include in the background model
 if hasattr(bpb3_config, 'max_ranking_for_background'):
   max_ranking_for_background = bpb3_config.max_ranking_for_background
 else:
   max_ranking_for_background=max_ranking_for_foreground

 if hasattr(bpb3_config,'max_ranking_for_significant_blocks'):
    max_ranking_for_significant_blocks=bpb3_config.max_ranking_for_significant_blocks
 else:
    max_ranking_for_significant_blocks=10

 resampling_iterations = bpb3_config.resampling_iterations
 gene_samples = bpb3_config.gene_samples
 background_shuffling_iterations = bpb3_config.background_shuffling_iterations
 background_command_file =  bpb3_config.background_command_file

 # 7. mutation signature definition for generating background mutations, set to None to disable
 mutation_signature_file = bpb3_config.mutation_signature_file  # "../../data/signature_7.tsv"
 # set to True to use tumor mutations from patients as background mutations, False otherwise.
 # This and the mutation signature are mutually exclusive.
 # In both are disabled, uniform mutations will be generated
 tumor_mutations_as_background =  bpb3_config.tumor_mutations_as_background

 # Results selection options for significant affinity changes
 # P value for ranksum test in the output
 affinity_change_p_value =  bpb3_config.affinity_change_p_value
 if hasattr(bpb3_config, 'pval_correction4affinity'):
      pval_correction4affinity=  bpb3_config.pval_correction4affinity
      if pval_correction4affinity:
         logging.info('Assign p-value correction in affinity change test')
      else:
         logging.info('No p-value correction in affinity change test')
 else:
      logging.info('Default False p-value correction in affinity change test')
      pval_correction4affinity=False

 logging.info(" %g "+ project_name + " patient expression datasets, %g normal expression datasets",len(donor_count_files),len(normal_count_files))


 #check step running conditions
 if hasattr(bpb3_config, 'runDiffExp'):
    runDiffExp=bpb3_config.runDiffExp
 else:
    runDiffExp=True

 if hasattr(bpb3_config, 'runRegionsFile'):
    runRegionsFile=bpb3_config.runRegionsFile
 else:
    runRegionsFile=True

 if hasattr(bpb3_config, 'runMussd'):
    runMussd= bpb3_config.runMussd
 else:
    runMussd=True

 if hasattr(bpb3_config, 'runHighMutBlocks'):
    runHighMutBlocks=bpb3_config.runHighMutBlocks
 else:
    runHighMutBlocks=True

 if hasattr(bpb3_config, 'runBayesPiBAR'):
    runBayesPiBAR=bpb3_config.runBayesPiBAR
 else:
    runBayesPiBAR=True

 if hasattr(bpb3_config, 'runBackgroundModel'):
    runBackgroundModel= bpb3_config.runBackgroundModel
 else:
    runBackgroundModel=True

 if hasattr(bpb3_config, 'runBlockSignificance'):
    runBlockSignificance=bpb3_config.runBlockSignificance
 else:
    runBlockSignificance=True

 if hasattr(bpb3_config, 'filterTFbyGene'):
    filterTFbyGene=bpb3_config.filterTFbyGene
 else:
    filterTFbyGene=True

 if hasattr(bpb3_config, 'isPlot'):
    isPlot=bpb3_config.isPlot
 else:
    isPlot=True

 if hasattr(bpb3_config, 'runRemoveTempFile'):
   runRemoveTempFile=bpb3_config.runRemoveTempFile
 else:
   runRemoveTempFile=True

 if hasattr(bpb3_config, 'runCluster4PWM'):
    runCluster4PWM=bpb3_config.runCluster4PWM
 else:
    runCluster4PWM=False

 if runCluster4PWM:
   #0.  Generate tempary pwms based on clustered and uncertain clustered PWMs ############
   if hasattr(bpb3_config, 'in_clustered_pwms_path') :
      in_clustered_pwms_path=bpb3_config.in_clustered_pwms_path
   else:
      in_clustered_pwms_path=out_data_folder
   
   if hasattr(bpb3_config, 'in_clustered_pwms_folder_name'):
     in_clustered_pwms_folder_name=bpb3_config.in_clustered_pwms_folder_name
   else:
     in_clustered_pwms_folder_name='quality_assessed_out_seed5'

   if hasattr(bpb3_config, 'in_uncertain_pwms_folder_name'):
     in_uncertain_pwms_folder_name=bpb3_config.in_uncertain_pwms_folder_name
   else:
     in_uncertain_pwms_folder_name='uncertain_pwms'

   if hasattr(bpb3_config, 'out_tmp_pwms_folder_name'):
     out_tmp_pwms_folder_name=bpb3_config.out_tmp_pwms_folder_name
   else:
     out_tmp_pwms_folder_name='tmp_pwm'
  
 # folder with PWM models for TF binding affinity
 #pwm_folder=os.path.join(out_data_folder,out_tmp_pwms_folder_name)
 pwm_folder = bpb3_config.pwm_folder


 ##########Below does not need to be changed for run the full pipeline ###############
 ############################## START BayesPI-BAR3 ###################################
 logging.info("")
 logging.info("          ---------------------------------")
 logging.info("         |                                 |")
 logging.info("         |  BayesPI-BAR3 pipeline started  |")
 logging.info("         |        |        |        |      |")
 logging.info("         |        V        V        V      |")
 logging.info("")
 logging.info("If the pipeline is interrupted at any point, simply restart it to continue the computations")
 logging.info("Pipeline parameters:")
 logging.info("\tDifferential gene expression: value is %s%s%sRPKM, P value for KS test: %g, minimum %sRPKM fold "
              "change: %s", "z-score transformed " * differential_expression_z_score_transform,
               "log-transformed " * differential_expression_log_transform,
               "quantile normalized " * differential_expression_quantile_normalization, differential_expression_p,
               "quantile normalized " * differential_expression_quantile_normalization,
               differential_expression_min_fold_change)

 logging.info("\tRegions of interest: [%d bp upstream --> TSS --> %d bp downstream]", upstream_size, downstream_size)
 logging.info("\tPatient data folder: %s", patient_data_folder)
 logging.info("\tControl gene expression folder: %s", normal_rnaseq_folder)
 logging.info("\tMuSSD options: hot region has at least %d mutations coming from at least %d different patients "
                 "within %d bp. %d bp added to both sides of each block.", min_block_size, min_patients_in_block,
                 cluster_distance, block_flank)

 logging.info("\tBayesPI-BAR options: %d sequence shuffling iterations, %d chemical potentials: %s",
                 shuffling_iterations, len(chemical_potentials.split(" ")), chemical_potentials)

 logging.info("\tBackground options: %d sequence shuffling iterations, including PWMs ranked at most %d in the "
                 "foreground, taking %d random regions %d times (%d total)", background_shuffling_iterations,
                 max_ranking_for_background, gene_samples, resampling_iterations, gene_samples * resampling_iterations)
 if pval_correction4affinity:
   logging.info("\tOutput options: Bonferroni-corrected Wilcoxon rank-sum test P value < %g", affinity_change_p_value)
 else:
   logging.info("\tOutput options: Wilcoxon rank-sum test P value < %g", affinity_change_p_value)

 logging.info("\tParallelization options are stored in the file %s: %s", parallel_options_file,
                 " ".join(common.read_lines(parallel_options_file)))

 #0. clustered PWMs
 if runCluster4PWM:
   logging.info("Step 0: make folder for clustered PWMs. BayesPI-BAR calculation and the top affected TFs will be based on these clustered PWMs!")
   cmd='bpb3 make_cluster4pwm --in_pwm_path ' + in_clustered_pwms_path +  \
        ' --in_clustered_pwm_folder ' +  in_clustered_pwms_folder_name + \
        ' --out_file_folder ' + out_tmp_pwms_folder_name + \
        ' --in_uncertain_pwm_folder ' + in_uncertain_pwms_folder_name

   if not os.path.exists(os.path.join(in_clustered_pwms_path,out_tmp_pwms_folder_name)):
     result_code=os.system(cmd)
     if result_code != 0:
        print("bpb3 make_cluster4pwm failed. See message in the console to find what went wrong.")
        exit(1)
     logging.info('Export at %s', os.path.join(in_clustered_pwms_path, out_tmp_pwms_folder_name))
   else:
     logging.info('Skip copy clustered PWMs because it exists in %s', os.path.join(in_clustered_pwms_path, out_tmp_pwms_folder_name))
   function2make_plots='make_plots_cluster4pwm'
   function2filter_results_by_gene_expression='filter_results_by_gene_expression_cluster4pwm'
 else:
   logging.info("Skip Step 0: make folder for clustered PWMs. BayesPI-BAR will use PWMs for the computation!")
   function2make_plots='make_plots'
   function2filter_results_by_gene_expression='filter_results_by_gene_expression'


 # 1. Differential expression
 logging.info("")
 if runDiffExp:
   logging.info("Step 1: Computing differentially expressed genes")

   cmd = "bpb3 differential_expression " + \
                                " --group_2_count_files " + " ".join(normal_count_files) + \
                                " --group_1_count_files " + " ".join(donor_count_files) + \
                                " --gene_lengths " + gene_length_file + \
                                " --output_file " + differentially_expressed_genes_file + \
                                " --output_group_1_rpkm " + group1_median_rpkm_file + \
                                " --output_group_2_rpkm " + group2_median_rpkm_file + \
                                " --p_value " + str(differential_expression_p) + \
                                " --test_method " + str(test_method) + \
                                " --quantile_normalization" * differential_expression_quantile_normalization + \
                                " --log_transform" * differential_expression_log_transform + \
                                " --z_score" * differential_expression_z_score_transform + \
                                " --output_all_values " * isExportAll +  \
                                (" --min_medianRPKM " + str(minimum_median_RPKM) if minimum_median_RPKM is not None else "") + \
                                (" --min_fold_change " + str(differential_expression_min_fold_change) if differential_expression_min_fold_change is not None else "")

   result_code= os.system(cmd)
   if result_code !=0:
       logging.info("Error in differential expresssion %s , ", cmd )
       exit(1)
  
   logging.info('Export at %s ', differentially_expressed_genes_file )
 else:
   logging.info('Skip step 1: differentially expressed genes')

 # 2. Make regions of interest near promoters of differentially expressed genes
 logging.info("")
 if runRegionsFile:
   logging.info("Step 2: Extracting regions near TSS of differentially expressed genes")

   cmd0="bpb3 gene_regions " +\
          " --gene_annotation " + gene_annotation_file +\
          " --transform_chrom" \
          " --output_file " + regions_file + \
          " --selected_genes " + differentially_expressed_genes_file + \
          " --upstream_size " + str(upstream_size) + \
          " --downstream_size " + str(downstream_size)
   result_code=os.system(cmd0)
   if result_code !=0:
       logging.info("Error in gene regions ")
       exit(1)
  
   logging.info('Export at %s ', regions_file )
 else:
   logging.info('Skip Step 2: generate region files  ')
   logging.info('	but use selected region file from input: %s', regions_file)

 # 3. Run MuSSD
 logging.info("")
 if runMussd:
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
                   (" --regions " + regions_file if regions_file is not None else "")
     #mussd_command_file=os.path.join(out_data_folder, 'run_mussd.sh')
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
   logging.info("Skip step 3: Mussd ")

 # 4. Select blocks that are mutated significantly more than average
 #significant_blocks_file = os.path.join(melanoma_folder, "significant_blocks.txt")
 logging.info("")
 if runHighMutBlocks:
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
   logging.info("Skip Step 4: selecting highly mutated blocks")

 significant_block_tags = frozenset(common.read_lines(significant_blocks_file))
 logging.info("%d blocks have significantly more mutations than expected: %s", len(significant_block_tags),
                 ", ".join(significant_block_tags))

 # 5. Compute foreground for each significant block
 logging.info("")
 if runBayesPiBAR:
   logging.info("Step 5: Running BayesPI-BAR on patient-specific mutation blocks")

   if not os.path.exists(foreground_folder) or start_from_integration :
     if not start_from_integration:
        os.mkdir(foreground_folder)

     for i, block_tag in enumerate(significant_block_tags):
        block_res_folder = os.path.join(foreground_folder, block_tag)
        ref_seq = os.path.join(mussd_result_folder, block_tag + "_ref.fasta")
        alt_seq = os.path.join(mussd_result_folder, block_tag + "_alt.fasta")
        if os.path.exists(block_res_folder):
            logging.info("The foreground output for block %s [%d of %d] already exists. BayesPI-BAR will now check its "
                         "consistency and recompute the broken or missing data ", block_tag, i + 1,
                         len(significant_block_tags))
        else:
            logging.info("BayesPI-BAR will now compute the delta-dbA scores and rankings for block %s [%d of %d]",
                         block_tag, i + 1, len(significant_block_tags))

        if start_from_integration:
          logging.info("BayesPI-BAR only starts from integration for %s  [%d of %d] ", block_tag, i+1, len(significant_block_tags) )
        result_code = os.system("bpb3 bayespi_bar " +
                                " --reference_sequences " + ref_seq +
                                " --alternate_sequences " + alt_seq +
                                " --pwm_folder " + pwm_folder +
                                " --chemical_potentials " + chemical_potentials +
                                " --iterations " + str(shuffling_iterations) +
                                " --result_folder " + block_res_folder +
                                " --max_rank " + str(max_ranking_for_foreground) +
                                " --reuse_output" + 
                                " --start_from_integration" * start_from_integration +
                                " @" + parallel_options_file)

#                                " --use_cores " + str(10) +
#                                " --max_nodes " + str(10) +
#                                " --use_slurm " + " --slurm_account " + 'nn4605k')

        if result_code !=0:
             logging.info("Error in bayespi bar ")
             exit(1)
        
        logging.info('Export at %s', block_res_folder )
   else:
     logging.info('BayesPI-BAR results already present, skipping step 5')
 else:
   logging.info("Skip Step 5: BayesPI-BAR ")

 # 6. Compute the background
 #background_command_file = os.path.join(melanoma_folder, "make_background.sh")
 logging.info("")
 if runBackgroundModel:
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
                            " --start_from_integration" * start_from_integration + 
                            " --max_rank " + str(max_ranking_for_background) +
                            ((" --mutation_signature " + mutation_signature_file)
                             if mutation_signature_file is not None else "") +
                            (" --background_mutations %s" %
                             os.path.join(mussd_result_folder, "mutations_in_regions.tsv")) *
                            tumor_mutations_as_background +
                            " @" + parallel_options_file)

#                            " --use_cores " + str(10) +
#                            " --max_nodes " + str(10) +
#                            " --use_slurm " + " --slurm_account " + "nn4605k")

#                            " @" + parallel_options_file)

   if result_code !=0:
       logging.info("Error in choose background parameters ")
       exit(1)

   logging.info('Export at %s ', background_command_file )

   #Failed at here !i
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
   logging.info("Skip Step 6: run background model ")

 # 7. Run significance tests
 logging.info("")
 if runBlockSignificance:
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
                                " --max_rank " + str(max_ranking_for_significant_blocks) +
                                " --p_value " + str(affinity_change_p_value) +
                                " --pval_correction " * pval_correction4affinity)
        if result_code !=0:
            logging.info("Error in affinity change significance test %s", block_tag)
            exit(1)
        
        logging.info("The table of significantly affected PWMs for block %s is stored in the file %s ", block_tag,
                     result_file)
 else:
   logging.info("Skip Step 7: block mutation significance test")

 # 8. Filter results by TF gene expression
 logging.info("")
 if filterTFbyGene:
   logging.info("Step 8: PWMs corresponding to the non-expressed TFs will now be removed from the results")

   for block_tag in significant_block_tags:
        block_folder = os.path.join(foreground_folder, block_tag)
        result_file = os.path.join(block_folder, "result.tsv")
        filtered_result_file = os.path.join(block_folder, "result_filtered.tsv")

        print("Filtering results by gene expression for block", block_tag, "...")
        result_code = os.system("bpb3 " + function2filter_results_by_gene_expression  + " "
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
 else:
   logging.info("Skip Step 8: filtering TF by gene expression ")
   for block_tag in significant_block_tags:
        block_folder = os.path.join(foreground_folder, block_tag)
        result_file = os.path.join(block_folder, "result.tsv")
        filtered_result_file = os.path.join(block_folder, "result_filtered.tsv")
        if os.path.isfile(result_file):
          print("Copy results by gene expression for block", block_tag, "...")
          result_code = os.system("\cp -r " + result_file +
                                " "  +  filtered_result_file )
          if result_code != 0:
            logging.info("Copying results for block %s failed. "
                         "See messages in the console to find what went wrong.", block_tag)
            exit(1)

 #9. plot all heatmaps
 logging.info("")
 if isPlot:
   logging.info("Step 9: Generate plots or heatmaps for the results")

   result_code=os.system('bpb3 ' + function2make_plots +' --foreground_folder ' + foreground_folder +
          ' --background_folder ' + background_folder + 
          ' --mussd_result_folder ' + mussd_result_folder )
   #bpb3 make_plots --foreground_folder '../../data/skin_cancer/out/foreground'  --background_folder '../../data/skin_cancer/out/pwm_for_background' --mussd_result_folder '../../data/skin_cancer/out/mussd_blocks'

   if result_code !=0:
       logging.info("Error in make plots ")
       exit(1)
   
   logging.info('Export at %s ', foreground_folder )
 else:
   logging.info("Skip step 9: plot heatmap")

 #10. clean temporary files
 logging.info("")
 if runRemoveTempFile:
   logging.info("Step 10: Remove temporary files in both background and foreground calculations. ")
   cmd= 'bpb3 clean_tmp --background_out_file ' + background_folder + \
                ' --foreground_out_file ' + foreground_folder  +  ' --chemical_potentials ' + chemical_potentials +  \
                ' --background_pwm_file ' + background_pwm_folder + ' --clean_both 2' 

   result_code=os.system(cmd )
   if result_code !=0:
     logging.info("Error in clean tmp ")
     logging.info(cmd )
     exit(1)
 else:
    logging.info("Skip step 10: remove temporary files")


#python clean_tmp.py --background_out_file '../../data/skin_cancer/out/background' \
#          --foreground_out_file '../../data/skin_cancer/out/foreground' \
#          --background_pwm_file '../../data/skin_cancer/out/pwm_for_background'


 logging.info("")
 logging.info("         |                                 |")
 logging.info("         |        |        |        |      |")
 logging.info("         |        V        V        V      |")
 logging.info("         |  BayesPI-BAR3 pipeline finished |")
 logging.info("          ---------------------------------")


if __name__ == "__main__":
  args=my_parser(argparse.ArgumentParser('run_pipeline.py')).parse_args()
  run(args)


