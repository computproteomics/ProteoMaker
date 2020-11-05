# PhosFake

This repository contains the scripts necessary for the generation of an in-silico bottom-up phosphoproteomics data set. 

All the parameters that are used to generate the data are listed in `Parameter.R`. The script that runs the entire pipeline is `RunAll.R`. Full batches of data set can be run via `FullBatchAnalysis.R`.

The pipeline is described in the figure `PhosFakeLayout.svg` and can be described as follows:

1) Generation of ground truth data at the proteoform level (`01_GenerateGroundTruth.R`).
2) Digestion of the proteoforms from the Ground Truth table (`02_Digestion.R`).
3) In silico MS run (`03_MSRun.R`).
4) Functions for data analysis from the peptide to proteins: `04_DataAnalysis.R`.
5) Statistical testing: `05_Statistics.R`
6) Benchmarking: `06_Benchmarks.R`
The folder `Benchmarking` contains the scripts used for proof of concept and benchmarking. 

### List of benchmarks

The following values and distribution are collected and will be used for comparing the results

#### Peptide level
- [X] Total number of peptides counting both modified and non-modified (fixed PTMs are not considered "modified")
- [X] Number of proteins (total of all proteins, independently whether shared or uniquely identified) 
- [X] Proportion of unique modified and non-modified peptides 
- [X] Total number of unique stripped sequences
- [X] Percentage missingness
- [X] AUC of ROC curve for correct differentially regulated features
- [X] TPR (found proportion of correct differentially regulated peptides, FDR < 0.05/0.01)
- [X] True FDR for estimated FDR of 0.01/0.05
- [X] Proportion of cleaved peptides (per number of miscleavages)
- [X] Retention time range (max-min)
- [X] Dynamic range (max-min on log2-scale)
- [X] Number of accepted PSMs (count scan numbers)
- [ ] TODO: Sum of squares residuals towards actual fold-changes (assigns mean fold-change for peptides with multiple assigned fold-changes )
- [ ] TODO: (only in silico) mean of std. dev. within replicates of peptides with actual fold-changes (on log-scale)
- [X] Distribution of quantitative values, how assymmetric (skewness) 


#### Protein level
- [X] Number of quantified protein groups
- [X] Proportion of quantified uniquely identified proteins
- [X] Percentage missingness
- [X] Distribution (summarized by its mean) of non-modified peptides per protein
- [X] AUC of ROC curve for correct differentially regulated features
- [X] True FDR for estimated FDR of 0.01/0.05
- [X] TPR (found proportion of correct differentially regulated proteins (FDR < 0.05/0.01)
- [X] Sum of squares residuals towards actual fold-changes
- [X] TODO: Dynamic range (max-min on log-scale)
- [X] TODO: (only in silico) mean of std. dev. within replicates (on log-scale)
- [X] Proportion of proteins with miscleaved peptides
- [X] Proportion of regulated proteins with wrong identified peptides (FDR < 0.01/0.05)
- [X] Distribution of quantitative values, how assymmetric (skewness) 

#### PTM level

 - [X] Number of quantitatively represented proteoforms in data set (including combined ones) (TODO: is that calculated from experimental data?)
 - [X] Number modified peptides
 - [X] Distribution of proteoforms per protein (summarized by mean?)
 - [X] Proportion of modified peptides with identical non-modified form
 - [X] AUC of ROC curve for correct differentially regulated modified peptides (after adjusting for protein amount)
 - [X] True FDR for estimated FDR of 0.01/0.05
 - [X] TPR (found proportion of correct differentially regulated proteins (FDR < 0.05/0.01)
 - [X] Proportion of wrongly regulated modified peptides (FDR 0.01/0.05, no adjustment)
 - [X] Proportion of modified peptides with quantified non-modified protein
 - [X] TODO: sum of square of residuals towards actual fold-changes (only modified peptides)
 

### Running full batches and benchmarking

The file `FullBatchAnalysis.R` allows running full batches without having to re-run the data sets which have been built with the same set of parameters. In addition, the pipeline is run hierarchically to avoid repetitive execution of identical down-stream analysis. This is done via creating hashes of the parameter configurations and writing intermediate and final results into respective tables.

Re-running the full batch with different assessment of the benchmarking metrics will avoid re-running the data set generation and analysis, and thus should be superfast.

Important remarks:

- You need to always define _all_ parameters in the beginning of this script
- Be aware that changing downstream parameters (ground truth, digestion) can immensely increase the number of possbile parameter settings
- Keep always the result files in the respective folder (`resultFilePath`) if you didn't change anything in the pipeline such as any of the methods in the sourced files
- For data set with only few quantified proteins (<100), benchmarking is discarding.


### Using PhosFake on the UCloud

Just mount the PhosFake folder into the application (RStudio, Jupyter, bash, ...) and work with a copy of `FullBatchAnalysis.R`. All new result files you generate will then be added and be available to the others.                   
