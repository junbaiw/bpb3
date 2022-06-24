---
title: 'BPB3: BayesPI-BAR version 3 package User Guide'
author: "Junbai Wang, Mingyi Yang, Magnar Bjoras, ..."
date: "`r format(Sys.time(), '%d %B, %Y')`"
output:
  html_document:
    df_print: paged
  word_document: default
---

<!--- R markdown set up --->
```{r setup, include=T, echo=FALSE, error=T, message=F}
knitr::opts_chunk$set(echo = TRUE)
library(rmarkdown)
library(stargazer)
library(reticulate)
reticulate::use_python("/Users/mingyi/opt/anaconda3/envs/snakemake/bin/python")
knitr::opts_chunk$set(engine.path = list(python = '/Users/mingyi/opt/anaconda3/bin/python'))
```

# Install package BPB3
**Down load package and unzip**
```{bash,  echo=T}
# tar -zxf bpb3.tgz
# cd bpb3
# move path to bpb3 fold, e.g: /cluster/projects/acount/user/bpb3_package/bpb3
# python setup.py install
```

**Readme in bpb3**

*BayesPI-BAR3*  
Type "python setup.py install" to install bpb3  package on your local machine.

final demo and data in folder: 
final_demo, final_demo_data

**Reference papers:**

1. Junbai Wang, Kirill Batmanov, BayesPI-BAR: a new biophysical model for characterization of regulatory sequence variations, Nucleic Acids Research, Volume 43, Issue 21, 2 December 2015, Page e147. [Link of Paper 1](https://academic.oup.com/nar/article/43/21/e147/2468103)
2. Batmanov, K., Wang, W., Bjørås, M. et al. Integrative whole-genome sequence analysis reveals roles of regulatory mutations in BCL6 and BCL2 in follicular lymphoma. Sci Rep 7, 7040 (2017). [Link of Paper 2](https://www.nature.com/articles/s41598-017-07226-4)
3. Batmanov K, Delabie J, Wang J. BayesPI-BAR2: A New Python Package for Predicting Functional Non-coding Mutations in Cancer Patient Cohorts. Front Genet. 2019. [Link of Paper 3](https://www.frontiersin.org/articles/10.3389/fgene.2019.00282/full)
4. Farooq A, Trøen G, Delabie J, Wang J. Integrating whole genome sequencing, methylation, gene expression, topological associated domain information in regulatory mutation prediction: A study of follicular lymphoma. Comput Struct Biotechnol J. 2022 Mar 23;20:1726-1742. [Link of Paper 4](https://www.sciencedirect.com/science/article/pii/S2001037022000976)   
5. Reference to our current paper where clustered PWMs are used ………


**Notes:**

1. Demo1 and demo2 were used in papers 1 and 3.
2. Demo3, demo4, demo5, and demo6 were used in papers 2 and 3.
3. Demo7 was similar to the calculation in paper 4, but the 14 FL samples were replaced by an independent FL cohort from papers 2 and 3.
4. Under final_demo_data, the pwm folder contains 1772 human pwms that were used in papers 1, 2, 3, and 4.
5. Applications of clustered PWMs in demos 2, 4, and 5 are described in current paper 5.  

# Demo1: Ranking the effect of sequence variation on TF binding by using 67 known regulatory SNPs

**demo1_67snps**

## Aim and pipeline flow
**Aim**
The goal of this demo is to predict putative target transcription factors (TFs) whose binding affinity is affected by the mutation at DNA sequence. In other words, to rank the significance of TF binding affinity changes due to the mutation (e.g., nucleotide polymorphisms -SNPs) in regulatory regions. More description of 67 regulatory SNPs that affect TF binding in this demo is available in paper 1.

**Pipeline**
![Figure 1. A computational pipeline of the Bayesian method for protein–DNA interaction with binding affinity ranking (BayesPI-BAR). The method is proposed for quantification of the effect of sequence variations on the binding affinity with transcription factors (TFs)](/Users/mingyi/course/R_markdown/data_bpb3/demo1_pipeline.png)

## Principle
Here, a published computational pipeline, a Bayesian method for protein–DNA interaction with binding affinity ranking (BayesPI-BAR), is used to quantify the effect of sequence variations on protein binding. The method includes two new parameters (TF chemical potentials or protein concentrations and differential TF binding affinity) that are not used in the previous application. This method is verified on 67 known human regulatory SNPs, of which 51 (76%) have predicted true TFs ranked in the top 10.

A biophysical modeling of protein-DNA interaction is used to estimate the binding probability between protein and DNA sequences, which is calculated by the Fermi–Dirac formula ([Djordjevic M et al, 2003](https://genome.cshlp.org/content/13/11/2381.full.pdf)): 
$$P(S)=1/(1+ e^{( E(S)−μ)/k*T })$$ 

where S represents the DNA sequence to be bound by a protein (i.e. 200 bp DNA sequences centered at ChIP-seq called peaks), E(S) is the protein binding free energy in sequence S or the correlated score in a *position-specific weight matrix* (PWM) of a TF; *mu* (μ) is the *chemical potential* (or protein concentration); k is a gas constant and T is the absolute temperature.

<!--- question: how to design /calculate the score for S ?, what is the unit? --->

**Estimation of chemical potential (*mu*):** The TF *chemical potential mu* at various conditions or cell lines can be estimated by fitting a known PWM to an in vivo ChIP-seq experiment, using Bayesian nonlinear regression model in package "*BayesPI2+*". Nevertheless, the *mu* can also be manually adjusted for a given PWM, such as by setting at range from 0 to -23 if the CHIP-deq experiments are missing. 

**Differential binding affinity (dbA)**:  This score is used to distinguish between the *direct binding* and the indirect protein–DNA interactions. The score of dbA is the difference between TF binding probability P(S) at the target DNA sequence and the binding probability at randomly mutated sequences. An expected P-value of dbA is used to remove indirect TF–DNA interactions after computing the TF binding affinity at a DNA sequence based on a set of collected TF PWMs.

**Shifted differential binding affinity**: After removing the putative indirect TF-DNA interaction from a set of collected PWMs,  the difference in TF binding affinity between the reference(Si, reference) and the mutated (Si, mutated) DNA sequence is quantified by an additional *shifted differential binding affinity* (δdbA), where δdbAi = dbA(Si,reference) − dbA(Si,mutate). 

**Ranking the effect of sequence variation (e.g., SNPs) on TF binding:** Usually, an exact *mu* of a given TF under a specific condition is usually unknown. The δdbA will be calculated by manually setting a range of *mu*.  To resolve this challenge, a principal component analysis (PCA) method was developed to combine the predictions from multiple chemical potentials (i.e. μ equals 0, −10, −13, −15, −18 and −20 in the present study). Then, it automatically ranks TFs based on their binding affinity changes caused by a regulatory sequence variation. The hypothesis is that the TFs with the strongest affected binding affinity change, due to a sequence variation, will be the one with the largest absolute value of δdbA(μ). Thus, BayesPI-BAR use the distance (or projected) scores for the final TF ranking.

More description of BayesPI2+ and Differential binding affinity (dba) are described in  [Junbai et al, 2016](https://pubmed.ncbi.nlm.nih.gov/26099425/).


## Run Job 1. Ranking the effect of sequence variation on TF binding
The example job is run a cluster computer (SAGA) provided by Norwegian Research Infrastructure Services(NRIS)

```
sbatch job_67snps.sm
```

```{bash,  echo=T, eval = FALSE}
#!/bin/csh

#Job name:
#SBATCH --job-name=demo1-bpb3
#
# Project:
#SBATCH --account=nn4605k
#
# Wall clock limit:
#SBATCH --time=12:00:00
#av
# Max memory usage:
#SBATCH --mem-per-cpu=20G --partition=bigmem
#
# Number of cores:
#SBATCH --cpus-per-task=12

#run 67 snps demo with PWMs
python bpb3_on_67snps_may22.py
```


**Data in**  
1. final_demo_data/demo1_67snps/in/snps_alt.fasta  
2. final_demo_data/demo1_67snps/in/snps_ref.fasta   
3. reference PWMs in original format

The set of 67 regulatory mutations (67 SNPs) were collected from HGMD database and other source, including 20 regulatory SNPs know to reduce or enhance the TF binding.

The position-specific weight matrices (PWMs) database were downloaded from earlier publications and resource such as JASPAR, containing over 2000 motifs representing over 600 unique human TFs. 

head of snps_alt.fasta,  the tile (annotation) for each SNP is in format "gene_TF|GAIN" or "gene_TF|LOSS".
```
>AGTRL1_SP1|GAIN
GTGCTCTCCGCCCTCCTGTTCTCACCCCCTCCCATCCAATCTAAATGGGACACATTATGCA
>ITGA2_SP1|LOSS
CCGGTGTTTGCGGAATCAGGAGGGGCGGGCCGGGGCGGGCCCTCGGCGCTGCAGGAGCTGC
>NFKBIL1_USF1|LOSS
```

head of snps_ref.fasta
```
>AGTRL1_SP1|GAIN
GTGCTCTCCGCCCTCCTGTTCTCACCCCCTTCCATCCAATCTAAATGGGACACATTATGCA
>ITGA2_SP1|LOSS
CCGGTGTTTGCGGAATCAGGAGGGGCGGGCTGGGGCGGGCCCTCGGCGCTGCAGGAGCTGC
>NFKBIL1_USF1|LOSS
GGCTGGAGGAAATGGCGCAAGCAGAGACGCAGGTGGAGGACGGAAGTGAACTGTGAGGGGC
```

<!-- ![Figure alignment of SNPS in SP1 with genome reference sequence](data_bpb3/alignment_SP1.png) -->
<!-- insert figurey by knit -->
```{r, screenshot, fig.align="left", out.width="60%", fig.cap="Alignment of SNPS in SP1 with genome reference sequence. \\label{sequence alignment}", echo=F}

knitr::include_graphics("data_bpb3/alignment_SP1.png")
```



**Output**

1. final_demo_data/demo1_67snps/out/demo1_67snps_bpb3.log
2. final_demo/demo1_67snps/out/foreground/test_seq/rankings/*.ranking.tsv      

Example of the output in .log file:
```
# Wed, 08 Jun 2022 20:05:10 INFO     Pipeline parameters:
# Wed, 08 Jun 2022 20:05:10 INFO     BayesPI-BAR options: 10000 sequence shuffling iterations, 6 chemical potentials: none -10 -13 -15 -18 -20
# Wed, 08 Jun 2022 20:33:53 INFO     Export at /cluster/projects/nn4605k/mingyi/bpb3_package/bpb3/final_demo_data/demo1_67snps/out/foreground/test_seq
# Wed, 08 Jun 2022 20:33:53 INFO     
# Wed, 08 Jun 2022 20:33:53 INFO     Remove temporary files in background or foreground calculations. 
```

Example of the output in TFs ranking:

* out file: ../../final_demo_data/demo1_67snps/out/foreground/test_seq/rankings/AGTRL1_SP1_GAIN_ranking.tsv

*Table* head of above output file, the TFs ranking for the sequence variation in gene AGTRL1.  
Please notice that 4 of the top 10 TFs are identified in the true target of transcription factor SP1:
```{r, echo=F}
d <- read.table("data_bpb3/AGTRL1_SP1_GAIN_ranking.tsv", header = T)
d <- head(d, 10)
knitr::kable(list(d))
```

**Main process**
The main processes are stored in Slurm file, example:
```
*Collecting dbA values...  
*Calculating delta-dbA values...  
*Integrating delta-dbA values across chemical potentials using mean_ddba ...  
*Integrating delta-dbA values across chemical potentials using pca ...  
*Computing PWM rankings for each variant..
```

**CPU requirement**  
  
* MaxRSS: 81 MB   
* Running time: 43 min    


## Run Job 2. Check predition accuracy

if you want to check the top N (e.g., 10) predicted TFs, which contain the true TF targets of the 67 regulatory SNPs, then please run the following script in terminal:
```{bash, eval=F}
bash check_accuracy4snps.sh 10
```

bash script: check_accuracy4snps.sh
```
#!/bin/bash
topN=$1
bpb3 check_accuracy --in_bpb3_ranking_path ../../final_demo_data/demo1_67snps/out/foreground/test_seq/rankings \
	--in_topRank_cutoff ${topN}
```

**output**
```
# Total : 67
# Considering top  10  TFs , change direction: False
# fredrikssonMuts
# 4
# 2  of  4  are correct ( 50.0 % )
# Mean rank:  1.0  Median rank: 1.0
# hgmdMuts
# 29
# 21  of  29  are correct ( 72.41379310344827 % )
# Mean rank:  3.0952380952380953  Median rank: 3.0
# andersonEpsteinMuts
# 34
# 28  of  34  are correct ( 82.35294117647058 % )
# Mean rank:  2.4642857142857144  Median rank: 1.5
# 
# 
# Total:  51  of  67  are correct ( 76.11940298507463 % )
# All: mean rank= 2.6666666666666665  median rank= 2.0
# Export at:  ../../final_demo_data/demo1_67snps/out/foreground/test_seq/rankings/../ranks.tsv
```

**output file**  
In total 67 SNPs, 51 have the true targeted TFs in the top 10 predictions. Among them, 21 SNPs with the true targeted TFs ranked in top 1.

```{r, echo= F, message=F}
d <- read.table("data_bpb3/ranks.tsv", header = T)
library(dplyr)
d <- d %>% arrange(bestRank)
#top10_count = nrow(d[d$bestRank %in% c(1:9), ]) # count of right prediction in the top 10 TFs, total 51 are correct.
#nrow(d[!is.na(d$bestRank), ])                  # alternative method for top10_count
#top1_count = nrow(d[d$bestRank %in% c(1), ])
knitr::kable(list(d[1:10, ]))
```


## Python Code

* bpb3_on_67snps_may22.py

```{python, echo = T, eval=F}
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
```

# End
-----------------------------------

