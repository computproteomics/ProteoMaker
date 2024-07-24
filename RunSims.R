################################################################################
#                       Run PhosFake                                           #
################################################################################

#####################
## Install PhosFake libraries
#####################
# install_phosfake()

#####################
## Load main functions
#####################
# You need to set the path to where the PhosFake (and this file) package is located
# get path of this file
path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)
source("R/00_BatchRunFuncs.R")

#####################
## Load PhosFake libraries and source files
#####################
load_phosfake()


#####################
## Paths and directories
#####################
phosfake_config <- set_phosfake(fastaFilePath = "Proteomes", resultFilePath = "QC_DataAnalysis/QC_output", cores = 8, clusterType = "FORK", calcAllBenchmarks = T)

#####################
## Create default list of testing parameters
#####################
Param <- def_param()

# Overwrite the default values with the ones you want to test
# You can use multiple values for each parameter that then will be combined for
# all possible combinations in different simulated datasets
# Param$paramGroundTruth$FastaFile <- "fasta_full_yeast.fasta"
Param$paramGroundTruth$NumReps <- c(3:5)
# Param$paramGroundTruth$NumCond <- 2
# Param$paramProteoformAb$QuantNoise <- seq(0.1, 0.9, 0.5)
# Param$paramProteoformAb$DiffRegFrac <- c(0.1, 0.3, 0.5)
# Param$paramProteoformAb$DiffRegMax <- seq(0.5, 2, 0.5)
# Param$paramDigest$Enzyme <- "trypsin"
# Param$paramDigest$PropMissedCleavages <- 0.01
# Param$paramDigest$MaxNumMissedCleavages <- 4
# Param$paramDigest$PepMinLength <- 7
# Param$paramDigest$PepMaxLength <- 30
# Param$paramMSRun$PercDetectedPep <- seq(0.1, 0.5, 0.1)


#####################

#####################
## Run the simulations
#####################
allBs <- run_sims(Param, phosfake_config)

#####################
## Get the results of an individual simulation
#####################
res <- get_simulation(allBs[[1]]$Param, phosfake_config)

#####################
## Make matrix of benchmarks and save
#####################
benchmarks <- matrix_benchmarks(allBs, phosfake_config)
write.csv(benchmarks, file = paste0(phosfake_config$resultFilePath, "/allBenchmarks.csv"))

#####################
## Visualize the results
#####################
visualize_benchmarks(benchmarks)
