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
- [ ] Number of peptides counting both modified and non-modified
- [ ] Number of proteins (total of different proteins, shared count as multiples)
- [ ] Proportion of unique peptides
- [ ] Proportion of shared peptides
- [ ] Percentage missingness
- [ ] AUC of ROC curve for correct differentially regulated features
- [ ] TPR (found proportion of correct differentially regulated peptides, FDR < 0.05/0.01)
- [ ] True FDR for estimated FDR of 0.01/0.05
- [ ] Proportion of miscleaved peptides
- [ ] Proportion of cleaved peptides (per number of miscleavages)

#### Protein level
- [ ] Number of quantified protein groups
- [ ] Proportion of quantified uniquely identified proteins
- [ ] Percentage missingness
- [ ] Distribution (summarized by its mean) of non-modified peptides per protein
- [ ] AUC of ROC curve for correct differentially regulated features
- [ ] True FDR for estimated FDR of 0.01/0.05
- [ ] TPR (found proportion of correct differentially regulated proteins (FDR < 0.05/0.01)
- [ ] Sum of squares residuals towards actual fold-changes
- [ ] Proportion of proteins with miscleaved peptides
- [ ] Proportion of regulated proteins with wrong identified peptides (FDR < 0.01/0.05)
- [ ] Distribution of quantitative values, how assymmetric (skewness) 

#### PTM level

 - [ ] Number of different proteoforms
 - [ ] Distribution of proteoforms per protein (summarized by mean?)
 - [ ] Number of peptides with modified and non-modified form
 - [ ] AUC of ROC curve for correct differentially regulated modified peptides (after adjusting for protein amount)
 - [ ] Proportion of wrongly regulated modified peptides (FDR 0.01/0.05, no adjustment)
 - [ ] Proportion of identifiable proteoforms (1+number of different modified peptides) with respect to actual proteoforms
 - [ ] Proportion of modified peptides with quantified non-modified protein
 


