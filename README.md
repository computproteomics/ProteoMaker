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
- [ ] Number of peptides 
- [ ] Number of proteins
- [ ] Number of unique peptides
- [ ] Number of shared peptides
- [ ] Percentage missingness
- [ ] Number of shared peptides
- [ ] AUC of ROC curve for correct differentially regulated features
- [ ] Number of differentially regulated peptides (FDR < 0.05/0.01)
- [ ] Percentage of miscleaved peptides

#### Protein level
- [ ] Number of quantified protein groups
- [ ] Number of quantified uniquely identified proteins
- [ ] Distribution (summarized by its mean) of peptides per protein
- [ ] AUC of ROC curve for correct differentially regulated features
- [ ] True FDR for estimated FDR of 0.01/0.05
- [ ] Number of differentially regulated proteins (FDR < 0.05/0.01)
- [ ] Sum of squares residuals towards actual fold-changes
- [ ] Percentage of proteins with miscleaved peptides
- [ ] Percentage of regulated proteins with wrong identified peptides (FDR < 0.01/0.05)

#### PTM level

 - [ ] Number of different proteoforms
 - [ ] Distribution of proteoforms per protein (summarized by mean?)
 - [ ] Number of peptides with modified and unmodified form
 - [ ] Number of 
 


