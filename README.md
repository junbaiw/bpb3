---
[comment]: <> (title: 'BPB3: BayesPI-BAR version 3 package User Guide')
[comment]: <> (author: "Junbai Wang, Mingyi Yang, Magnar Bjoras, ...")
[comment]: <> (date: "`r format(Sys.time(), '%d %B, %Y')`" )
[comment]: <> (output: )
[comment]: <> (  html_document: )
[comment]: <> (    df_print: paged )
[comment]: <> (  word_document: default)
# BPB3: BayesPI-BAR version 3 User Guide
---

<!--- R markdown set up --->
<!-- ```{r setup, include=T, echo=FALSE, error=T, message=F}
knitr::opts_chunk$set(echo = TRUE)
library(rmarkdown)
library(stargazer)
library(reticulate)
reticulate::use_python("/Users/junbai/miniconda3/bin//python")
knitr::opts_chunk$set(engine.path = list(python = '/Users/junbai/miniconda3/bin//python'))
```
-->

# Install package BPB3    
## Download package and unzip

```{bash,  echo=T}
# tar -zxf bpb3.tgz
# cd bpb3
# python setup.py install
```

## Readme in bpb3

*BayesPI-BAR3*  
Type "python setup.py install" to install bpb3  package on your local machine.        
file folders: final_demo, final_demo_data


## Reference papers

1. Junbai Wang, Kirill Batmanov, BayesPI-BAR: a new biophysical model for characterization of regulatory sequence variations, Nucleic Acids Research, Volume 43, Issue 21, 2 December 2015, Page e147. [Link of Paper 1](https://academic.oup.com/nar/article/43/21/e147/2468103)

2. Batmanov, K., Wang, W., Bjørås, M. et al. Integrative whole-genome sequence analysis reveals roles of regulatory mutations in BCL6 and BCL2 in follicular lymphoma. Sci Rep 7, 7040 (2017). [Link of Paper 2](https://www.nature.com/articles/s41598-017-07226-4)

3. Batmanov K, Delabie J, Wang J. BayesPI-BAR2: A New Python Package for Predicting Functional Non-coding Mutations in Cancer Patient Cohorts. Front Genet. 2019. [Link of Paper 3](https://www.frontiersin.org/articles/10.3389/fgene.2019.00282/full)

4. Farooq A, Trøen G, Delabie J, Wang J. Integrating whole genome sequencing, methylation, gene expression, topological associated domain information in regulatory mutation prediction: A study of follicular lymphoma. Comput Struct Biotechnol J. 2022 Mar 23;20:1726-1742. [Link of Paper 4](https://www.sciencedirect.com/science/article/pii/S2001037022000976)   

5. Reference to our current paper where clustered PWMs are used ………


## Notes

1. Demo1 and demo2 were used in papers 1 and 3.
2. Demo3, demo4, demo5, and demo6 were used in papers 2 and 3.
3. Demo7 was similar to the calculation in paper 4, but the 14 FL samples were replaced by an independent FL cohort from papers 2 and 3.
4. Under final_demo_data, the pwm folder contains 1772 human pwms that were used in papers 1, 2, 3, and 4.
5. Applications of clustered PWMs in demos 2, 4, and 5 are described in current paper 5.  

[More information about demos in bpb3 package please refer to ->  https://bpb3.github.io/bpb3 ](https://bpb3.github.io/bpb3/)

[bpb3 package can be downloaded from -> https://github.com/junbaiw/bpb3 ](https://github.com/junbaiw/bpb3)
# End


