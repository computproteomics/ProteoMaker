# PhosFake: Simulating proteoforms in bottom-up proteomics

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4540737.svg)](https://doi.org/10.5281/zenodo.4540737)

Phosfake is a platform for the  generation of an in-silico bottom-up proteomics data set with a ground truth on the level of proteoforms. 

All the parameters that are used to generate the data are described in [man/Parameters.qmd](man/Parameters.html). The script that runs the entire pipeline is [`RunSims.R`](inst/cmd/RunSims.R). Alternatively, you can use the [Vignette](vignettes/Vignette.html). The simulations with multiple parameters can be set up and run. PhosFake also provides comparison of the results with the ground truth using [benchmarking metrics](man/Benchmarks.html) and visual comparison between simulated data sets.

The pipeline is described in the figure [`PhosFakeLayout.svg`](inst/img/PhosFakeLayout) and can be described as follows:

0) General functions to run the simulations: [`00_BachRunFuncs.R`](R/00_BatchRunFuncs.R)
1) Generation of ground truth data at the proteoform leve:l [`01_GenerateGroundTruth.R`](R/01_GenerateGroundTruth.R).
2) Digestion of the proteoforms from the ground truth: [`02_Digestion.R`](R/02_Digestion.R).
3) In silico MS run [`03_MSRun.R`](R/03_MSRun.R).
4) Functions for data analysis from the peptide to proteins: [`04_DataAnalysis.R`](R/04_DataAnalysis.R).
5) Statistical testing: [`05_Statistics.R`](R/05_Statistics.R).
6) Benchmarking: [`06_Benchmarks.R`](R/06_Benchmarks.R).

## Installation

Install the package from GitHub with `devtools`:
`devtools::install_github("veitveit/PhosFake")`

### Running full batches and benchmarking

Running the [Vignette](vignettes/Vignette.html) allows running full batches without having to re-run the data sets which have been built with the same set of parameters. In addition, the pipeline is run hierarchically to avoid repetitive execution of identical down-stream analysis. This is done via creating hashes of the parameter configurations and writing intermediate and final results into respective tables.

Re-running the full batch with different assessment of the benchmarking metrics will avoid re-running the data set generation and analysis, and thus should be superfast.

Important remarks:

- You need to always define _all_ parameters in the beginning of this script
- Be aware that changing downstream parameters (ground truth, digestion) can immensely increase the number of possible parameter settings
- Keep always the result files in the respective folder (`resultFilePath`) if you didn't change anything in the pipeline such as any of the methods in the sourced files. This will allow you to run the full batch without re-running the data set generation and analysis.
- For data set with only few quantified proteins (<100), benchmarking is discarding.

