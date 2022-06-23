import argparse
import sys
#author Junbai Wang September 2021

class Main(object):
  def __init__(self):
     parser = argparse.ArgumentParser(
              description= "BayesPI-BAR in Python3 - bpb3",
              usage=''' bpb3 <task> [<args>]

     Tasks available for using:
         differential_expression	Predict differentialy expressed genes (DEG) based on two group of samples.
         gene_regions		Extracts regions near transcription start sites of selected genes based on genCode gtf. 
         mussd			Mutation filtering based on the Space and Sample Distribution - MuSSD.
         highly_mutated_blocks	Find blocks with significantly more mutations than would be expected.
         bayespi_bar		BayesPI-BAR delta-dbA ranking computation for TF binding affinity affected by DNA mutation.
         choose_background_parameters	Selects parameters for mutation background computation.
         background_affinity_changes	Mutation background computation.
         affinity_change_significance_test	Significant test of TF binding affinity changes between foreground and background affinity changes.
         parallel		Run commands from the given file in parallel.
         make_cluster4pwm	Make input PWM files for bpb3 based on clustered PWMs.
         bpb3selectedPWM	The second level analysis of bpb3 by using the top PWMs in TF ranking after the first level analysis of bpb3 based on the clustered PWMs. 
         run_pipeline		Run full bpb3 pipeline (e.g., the first level analysis of bpb3 if clustered PWMs are used in the calculation).
	 clean_tmp		Clean temporary files from output folders.

     Tasks available for demo purpose:
         plot_result 		Generate heatmaps for selected mutation blocks. (demo) 
         filter_results_by_gene_expression_cluster4pwm	Filter those TF whose expression is too low in clustered PWMs. (demo)
         filter_results_by_gene_expression		Filter those TFs whose expression is too low (e.g., RPKM<0.03). (demo)
         make_plots_cluster4pwm		Make heatmap plots for all significant mutation blocks that affecting clustered PWMs. (demo)
         make_plots			Make heatmap plots for all significant mutation blocks that affecting PWMs. (demo)
         check_accuracy4cluster	Check accuracy for 67 SNPs that based on clustered PWMs. (demo)
         check_accuracy		Check accuracy for 67 SNPs that based on original PWMs. (demo)
         filterDEG4bpb3		Filter bpb3 exported differential expression gene list by rratios. (demo)
         preprocess_icgc_data	Preprocess of ICGC data such as a folder contains files donor_*, specimen, simple_somatic*, exp_seq.tsv et al. (demo)
     ''' )

     parser.add_argument('task', help ='Pipeline task to run')
     args= parser.parse_args(sys.argv[1:2])
     if not hasattr(self, args.task):
          print('*******Error: Unrecognized task ********')
          parser.print_help()
          exit(1)
     getattr(self,args.task)()

  def differential_expression(self):
      from .script.differential_expression import my_parser, run 
      parser = my_parser(argparse.ArgumentParser(prog='bpb3 differential_expression',
               description="This program determines which genes are differentially expressed "
                           "based on RNA-seq data for two groups of samples. RPKM values are "
                           "computed for each sample, optionally normalized, and Kolmogorov-Smirnov test/T-test "
                           "is then applied to them to determine significant difference between "
                           "distributions of values of the two groups. Order of optional normalizations: "
                           "1) quantile normalization, 2) log transform, 3) z-score transform. "))
      run(parser.parse_args(sys.argv[2:]))

  def preprocess_icgc_data(self):
      from .script.preprocess_icgc_data import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='bpb3 preprocess_icgc_data',
               description="This is a demo script for getting and preprocessing of ICGC data before runnding bpb3 "))
      run(parser.parse_args(sys.argv[2:]))

  def gene_regions(self):
      from .script.gene_regions import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='bpb3 gene_regions',
               description="This program extracts regions near transcription start site of"
                            "selected genes based on gencode gtf file."))
      run(parser.parse_args(sys.argv[2:]))

  def mussd(self):
      from .script.mussd import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='bpb3 mussd',
               description="Mutation filtering based on the Space and Sample Distribution - MuSSD."))
      run(parser.parse_args(sys.argv[2:]))
  
  def highly_mutated_blocks(self):
      from .script.highly_mutated_blocks import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='bpb3 highly_muated_blocks',
               description="This program finds which blocks produced by mussd.py have"
                           "significantly more mutations than would be expected if all mutations "
                           "were uniformly distributed across regions of interest."))
      run(parser.parse_args(sys.argv[2:]))

  def bayespi_bar(self):
      from .script.bayespi_bar import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='bpb3 bayespi_bar',
               description="BayesPI-BAR delta-dbA ranking computation for TFs binding affinity changes affected by DNA mutations. "
                           "This program computes rankings of given PWMs/clustered PWMs according to predicted effects of "
                           "given sequence variants." ,fromfile_prefix_chars="@" ))
      run(parser.parse_args(sys.argv[2:]))

  def choose_background_parameters(self):
      from .script.choose_background_parameters import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='bpb3 choose_background_paramseters',
               description="""This program selects parameters (block size, mutation distribution, set of PWM files) 
                             for background computation based on selected mutation blocks obtained by mussd.py. 
                            It will create a shell script with the command to create the background model with selected parameters.""" , fromfile_prefix_chars="@"))
      run(parser.parse_args(sys.argv[2:]))

  def background_affinity_changes(self):
      from .script.background_affinity_changes import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='bpb3 background_affinity_changes',
                description="This program is used to compute mutation backgrounds based on selected parameters. ",
                fromfile_prefix_chars="@"))
      run(parser.parse_args(sys.argv[2:]))

  def affinity_change_significance_test(self):
      from .script.affinity_change_significance_test import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='bpb3 affinity_change_significance_test',
               description="Significant test of TF binding affinity changes between the foreground and the background calculations"))
      run(parser.parse_args(sys.argv[2:]))

  def parallel(self):
      from .script.script_high.parallel import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='bpb3 parallel',
               description="This program runs commands from the given file in parallel."))
      run(parser.parse_args(sys.argv[2:]))

  def plot_result(self):
      from .script.script_high.plot_result import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='bpb3 plot_resukt',
               description= "This demoe script produces heatmaps like the one in the package paper. For general " 
                            "use, you may need to edit the layout, font sizes etc. so that everything is visible.")) 
      run(parser.parse_args(sys.argv[2:]))

  def make_plots(self):
      from .script.make_plots import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='bpb3 make_plots',
               description="This is a demo script to generate heatmaps for all significant mutation blocks that affect TF binding based on calculation of PWMs."))
      run(parser.parse_args(sys.argv[2:]))
  
  def make_plots_cluster4pwm(self):
      from .script.make_plots_cluster4pwm import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog= 'bpb3 make_plots_cluster4pwm',
               description="This is a demo script to generate heatmaps for all significant mutation blocks that affect TF binding based on calculation of clustered PWMs"))
      run(parser.parse_args(sys.argv[2:]))

  def filter_results_by_gene_expression(self):
      from .script.script_high.filter_results_by_gene_expression import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='bpb3 filter_results_by_gene_expression',
               description= "This program will filter those TFs whose expression is too low (e.g., RPKM<0.03) from "
                            "the results produced by affinity_change_significance_test.py. "
                            "Warning: it assumes a certain naming convention for PWM files. "
                            "If you add your own PWM files that do not conform to it it will "
                            "likely crash or filter them out. It is also incomplete and doesn't "
                            "recognize many TFs. This is the reason it is not a part of " 
                            "the package, but rather a demo" ))
      run(parser.parse_args(sys.argv[2:]))
   
  def filter_results_by_gene_expression_cluster4pwm(self):
      from .script.script_high.filter_results_by_gene_expression_cluster4pwm import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='bpb3 filter_result_by_gene_expression_cluster4pwm',
                 description="This function will filter TFs with too low expresssion (e.g., RPKM<0.03) from the results "
                   "producted by affinity_change_significance_test.py for clustered PWMs. Warning : this program is for demo purpose"))
      run(parser.parse_args(sys.argv[2:]))

  def clean_tmp(self):
      from .script.script_high.clean_tmp import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='bpb3 clean_tmp',
               description="Clean temporary files in output folders, default only clean temporary files in background output. "
                           "This function will remove intermediate result files that are not used anywhere so that it will reduce the result file size. "))
      run(parser.parse_args(sys.argv[2:]))

  def run_pipeline(self):
      from .script.run_pipeline import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='bpb3 run_pipeline',
               description="Run BayeaPI-BAR3 full pipeline with a proper bpb3 configure file. " 
               "If clustered PWMs are used in the calcuation then it is the first level analysis of bpb3. "
               "The second level analysis of bpb3, which evaluates the top N PWMs from the ranking test of the first analysis, can be achieved by function bpb3selectedPWM."))
      run(parser.parse_args(sys.argv[2:]))
 
  def make_cluster4pwm(self):
      from .script.script_high.other.make_cluster4pwm import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog='bpb3 make_cluster4pwm',
              description='This function makes input PWMs files for bpb3 by using clustered PWMs from abc4pwm.  '))
      run(parser.parse_args(sys.argv[2:]))

  def check_accuracy4cluster(self):
      from .script.script_high.other.check_accuracy4cluster import my_parser, run
      parser= my_parser(argparse.ArgumentParser(prog= 'bpb3  check_accuracy4cluster', 
              description="This is a demo for checking accuracy of predicted results in 67 SNPs based on clustered PWMs."))
      run(parser.parse_args(sys.argv[2:]))

  def check_accuracy(self):
      from .script.script_high.other.check_accuracy import my_parser, run
      parser= my_parser(argparse.ArgumentParser(prog= 'bpb3 check_accuracy',
              description = "This is a demo for checking accuracy of predictions for 67 SNPs based on original PWMs."))
      run(parser.parse_args(sys.argv[2:]))

  def bpb3selectedPWM(self):
      from .script.bpb3selectedPWM import my_parser2, run
      parser= my_parser2(argparse.ArgumentParser(prog= 'bpb3 bpb3selectedPWM',
              description="This is the second level analysis of bpb3, to evaluate PWMs of the selected top N clustered PWMs from "
                           " the first level analysis of bpb3 by using the clustered PWMs.",
              fromfile_prefix_chars="@"))
      run(parser.parse_args(sys.argv[2:]))

  def filterDEG4bpb3(self):
      from .script.filterDEG4bpb3 import my_parser, run
      parser = my_parser(argparse.ArgumentParser(prog= 'bpb3 filterDEG4bpb3',
               description='This is a demo script to filter differential expression gene (DEG) list exported by bpb3 by using rratio and minimum RPKM in each group'))
      run(parser.parse_args(sys.argv[2:]))

def main():
  Main()

if __name__== '__main__':
  Main()


    
