################################################################################
#                       Run ProteoMaker                                           #
################################################################################


#####################
## Paths and directories
#####################
proteomaker <- set_proteomaker(fastaFilePath = system.file("Proteomes", package = "ProteoMaker"),
                                resultFilePath = "SimulatedDataSets",
                                cores = 2, clusterType = "PSOCK",
                                runStatTests = T,
                                calcAllBenchmarks = T
                                )

#####################
## Create default list of testing parameters
#####################
Param <- def_param("tmp/parameters.yml")
# Param <- def_param()


# Overwrite the default values with the ones you want to test
# You can use multiple values for each parameter that then will be combined for
# all possible combinations in different simulated datasets
Param$paramGroundTruth$PathToFasta <- "fasta_full_human.fasta"


Param$paramGroundTruth$NumReps <- c(5)
# Param$paramGroundTruth$NumCond <- 5
Param$paramGroundTruth$PercExpressedProt <- 0.5

Param$paramGroundTruth$FracModProt <- 0
# Param$paramGroundTruth$PTMTypes <- list(mods=c("ph", "ox"))
# Param$paramGroundTruth$PTMTypesDistr <- list(m1 = c(ph=0.5, ox = 0.4))
# Param$paramGroundTruth$PTMTypesMass <- list(m1 = c(ph=79.966331, ox = 15.994915))
# Param$paramGroundTruth$ModifiableResidues <- list(m1 = list(ph=c("S", "T", "Y"), ox = c("M", "K")))
# Param$paramGroundTruth$ModifiableResiduesDistr <- list(m1 = list(ph=c(0.86, 0.13, 0.01), ox = c(0.5, 0.5)))
# Param$paramGroundTruth$PTMTypes <- list(mods1 = c("ph", "ac"))
# Param$paramGroundTruth$PTMTypesMass <- list(m1 = c(ph=79.966331, ac=42.010565))
# Param$paramGroundTruth$PTMTypesDistr <- list(m1 = list(ph=0.5, ac=0.5), m2 = list(ph=0.1, ac=0.9))
# Param$paramGroundTruth$ModifiableResidues <- list(m1 = list(ph=c("S", "T", "Y"), ac=c("K", "R")))
# Param$paramGroundTruth$ModifiableResiduesDistr <- list(m1 = list(ph=c(0.86,0.13, 0.01), ac=c(0.5, 0.5)))

# Param$paramGroundTruth$PTMTypes <- list(mods=c("ph", "ox"))
# Param$paramGroundTruth$PTMTypesDistr <- list(m1 = c(ph=0.5, ox = 0.04))
# Param$paramGroundTruth$PTMTypesMass <- list(m1 = c(ph=79.966331, ox = 15.994915))
# Param$paramGroundTruth$ModifiableResidues <- list(m1 = list(ph=c("S", "T", "Y"), ox = c("M")))
# Param$paramGroundTruth$ModifiableResiduesDistr <- list(m1 = list(ph=c(0.86, 0.13, 0.01), ox = c(1)))
# Param$paramGroundTruth$PTMMultipleLambda <- c(0.5)
# Param$paramDataAnalysis$MinUniquePep <- 100

# Param$paramProteoformAb$QuantNoise <- c(0.5)
#Param$paramProteoformAb$DiffRegFrac <- c(0.1, 0.3, 0.5)
 Param$paramProteoformAb$DiffRegMax <- c(5)
# Param$paramDigest$Enzyme <- "trypsin"
# Param$paramDigest$PropMissedCleavages <- 0.01
# Param$paramDigest$MaxNumMissedCleavages <- 4
# Param$paramDigest$PepMinLength <- 7
# Param$paramDigest$PepMaxLength <- 30

# Param$paramDigest$EnrichmentEfficiency <- 0.8
Param$paramDigest$EnrichPTM <- NA
# Param$paramDigest$ModificationLoss <- 0.5
Param$paramMSRun$PercDetectedVal <- 0.5
Param$paramMSRun$PercDetectability <- 0.25
# Param$paramMSRun$WrongLocalizations <- 0.1
Param$paramMSRun$WrongIDs <- 0 #seq(0,0.3,0.01)

#####################
## Read parameters from yaml file
#####################r
# Alternatively, you can set and read the parameters from a yaml file
# See the example file in the main folder
# Param <- def_param_from_yaml("path/to/your/param.yaml")


#####################

#####################
## Run the simulations
#####################
allBs <- run_sims(Param, Config = proteomaker)

#####################
## Get the results of an individual simulation
#####################
res <- get_simulation(allBs[[1]]$Param, Config = proteomaker)

#####################
## Make matrix of benchmarks and save
#####################
benchmarks <- matrix_benchmarks(allBs, )
write.csv(benchmarks, file = paste0(proteomaker_conf$resultFilePath, "/allBenchmarks.csv"))

#####################
## Visualize the results
#####################
# visualize the benchmarks and parameters of the second simulation
visualize_benchmarks(benchmarks,2)


