################################################################################
#                       Run PhosFake                                           #
################################################################################
library(PhosFake)

#####################
## Paths and directories
#####################
phosfake_config <- set_phosfake(fastaFilePath = "Proteomes", resultFilePath = "SimulatedDataSets", cores = 8, clusterType = "FORK", calcAllBenchmarks = T)

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
Param$paramMSRun$DetectabilityThreshold <- 0.5
Param$paramMSRun$PercDetectedVal <- 0.75

#####################
## Read parameters from yaml file
#####################
# Alternatively, you can set and read the parameters from a yaml file
# See the example file in the main folder
# Param <- def_param_from_yaml("path/to/your/param.yaml"


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
# visualize the benchmarks and parameters of the second simulation
visualize_benchmarks(benchmarks, 2)

