# PhosFake

This repository contains the scripts necessary for the generation of an in-silico bottom-up phosphoproteomics data set. 

All the parameters that are used to generate the data are listed in `Parameter.R`. The script that runs the entire pipeline is `RunAll.R`.

The pipeline is described in the figure `PhosFakeLayout.svg` and can be described as follows:

1) Generation of ground truth data at the proteoform level (`01_GenerateGroundTruth.R`).
2) Digestion of the proteoforms from the Ground Truth table (`02_Digestion.R`).
3) In silico MS run (`03_MSRun.R`).
4) Functions for data analysis from the peptide to proteins: `04_DataAnalysis.R`.
5) Statistical testing
6) (TODO) Benchmarking
The folder `Benchmarking` contains the scripts used for proof of concept and benchmarking. 

### List of benchmarks

The following values and distribution are collected and will be used for comparing the results

#### Peptide level
- [X] Number of peptides counting both modified and non-modified
- [X] Number of proteins (total of different proteins, shared count as multiples)
- [X] Proportion of unique peptides
- [X] Proportion of shared peptides
- [X] Percentage missingness
- [X] AUC of ROC curve for correct differentially regulated features
- [X] TPR (found proportion of correct differentially regulated peptides, FDR < 0.05/0.01)
- [X] True FDR for estimated FDR of 0.01/0.05
- [X] Proportion of miscleaved peptides
- [X] Proportion of cleaved peptides (per number of miscleavages)

#### Protein level
- [X] Number of quantified protein groups
- [X] Proportion of quantified uniquely identified proteins
- [X] Percentage missingness
- [X] Distribution (summarized by its mean) of non-modified peptides per protein
- [X] AUC of ROC curve for correct differentially regulated features
- [X] True FDR for estimated FDR of 0.01/0.05
- [X] TPR (found proportion of correct differentially regulated proteins (FDR < 0.05/0.01)
- [X] Sum of squares residuals towards actual fold-changes
- [X] Proportion of proteins with miscleaved peptides
- [X] Proportion of regulated proteins with wrong identified peptides (FDR < 0.01/0.05)
- [X] Distribution of quantitative values, how assymmetric (skewness) 

#### PTM level

 - [X] Number of quantitatvely represented proteoforms in data set (including combined ones)
 - [X] Number modified peptides
 - [X] Distribution of proteoforms per protein (summarized by mean?)
 - [X] Proportion of modified peptides with identical non-modified form
 - [X] AUC of ROC curve for correct differentially regulated modified peptides (after adjusting for protein amount)
 - [X] True FDR for estimated FDR of 0.01/0.05
 - [X] TPR (found proportion of correct differentially regulated proteins (FDR < 0.05/0.01)
 - [X] Proportion of wrongly regulated modified peptides (FDR 0.01/0.05, no adjustment)
 - [X] Proportion of modified peptides with quantified non-modified protein
 


