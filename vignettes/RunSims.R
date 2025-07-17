## -----------------------------------------------------------------------------
# Install ProteoMaker
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = "https://cran.r-project.org")
}
devtools::install_github("computproteomics/ProteoMaker")

## -----------------------------------------------------------------------------
# Load ProteoMaker
library(ProteoMaker)

## -----------------------------------------------------------------------------
# Configure the paths and settings
proteomaker_config <- set_proteomaker(
  # Default path to the Proteomes folder in the ProteoMaker package
  fastaFilePath =  system.file("Proteomes", package = "ProteoMaker"),
  # Default path to the temporary directory
  resultFilePath = paste0(tempdir(), "/SimulatedDataSets"),
  # Increase to the number of cores you want to use.                                                                     
  # Large numbers can lead to considerably higher RAM usage
  cores = 1, 
  # Should work on all operation systems. If you have problems, try "FORK" instead
  clusterType = "PSOCK",
  # Run statistical tests on the simulated data to compare ground truth with simulated values
  runStatTests = TRUE,
  # Calculate all benchmarks for the simulated data
  calcAllBenchmarks = TRUE               
)

## -----------------------------------------------------------------------------
# Generate default parameters
Param <- def_param()


# Example of overwriting default values
Param$paramGroundTruth$PathToFasta <- "fasta_example.fasta"
Param$paramGroundTruth$PercExpressedProt <- 1.0
# Param$paramGroundTruth$NumReps <- c(3)
# Param$paramGroundTruth$NumCond <- 5
# Param$paramGroundTruth$FracModProt <- 0.5
# Param$paramGroundTruth$PTMTypes <- "ph"
# Param$paramGroundTruth$PTMTypesMass <- c(79.966331)
# Param$paramGroundTruth$PTMTypesDist <- c(1)
# Param$paramGroundTruth$PTMMultipleLambda <- c(0.1)
# Param$paramGroundTruth$ModifiableResidues <- list(c("S", "T", "Y"))
# Param$paramGroundTruth$ModifiableResiduesDistr <- list(c(0.86,0.13, 0.01))
Param$paramGroundTruth$NumReps <- c(3:5)
# Param$paramProteoformAb$QuantNoise <- seq(0.1, 0.9, 0.5)
# Param$paramProteoformAb$DiffRegFrac <- c(0.1, 0.3, 0.5)
# Param$paramProteoformAb$DiffRegMax <- seq(0.5, 2, 0.5)
# Param$paramDigest$Enzyme <- "trypsin"
# Param$paramDigest$PropMissedCleavages <- 0.01
# Param$paramDigest$MaxNumMissedCleavages <- 4
# Param$paramDigest$PepMinLength <- 7
# Param$paramDigest$PepMaxLength <- 30
# Param$paramMSRun$PercDetectedPep <- seq(0.1, 0.5, 0.1)
# Param$paramDataAnalysis$MinUniquePep <- 100


## -----------------------------------------------------------------------------
# Run the simulations
allBs <- run_sims(Param, proteomaker_config)


## -----------------------------------------------------------------------------
# Retrieve results}
res <- get_simulation(allBs[[1]]$Param, proteomaker_config)


## -----------------------------------------------------------------------------
# Generate the benchmark matrix
benchmarks <- matrix_benchmarks(allBs, proteomaker_config)
write.csv(benchmarks, file = paste0(proteomaker_config$resultFilePath, "/allBenchmarks.csv"))


