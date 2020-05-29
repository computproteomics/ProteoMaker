################################################################################
#                       QC batch test of the parameters                        #
################################################################################

library(purrr)
library(crayon)
library(parallel)
library(digest)

## only once for installation
# myPaths <- .libPaths()   # get the paths
# .libPaths("/home/rstudio/R")
# library("BiocManager")
# install(c("protr","preprocessCore","matrixStats","extraDistr","fdrtool","qvalue","limma","moments"))

#### YOU NEED TO BE IN THE MAIN FOLDER OF THE PHOSFAKE SCRIPTS

####################### Paths and directories
#####################
#Working directory should be PhosFake's main directory.
# File paths
Param <- NULL
fastaFilePath <- "QC_DataAnalysis"
resultFilePath <- "QC_DataAnalysis/QC_output"
# For parallel computing
cores <- 20
clusterType <- "FORK"
calcAllBenchmarks <- T

try(  dir.create(resultFilePath))


repportName <- paste0(resultFilePath,"/QC_Report_BatchAnalysis_", gsub(":|[[:space:]]", "_", format(Sys.time(), "%Y%b%d%X")))

file.create(paste0(repportName,".txt"))
#####################

#####################
## Load sources
#####################

source("Parameter.R")
source("01_GenerateGroundTruth.R")
source("02_Digestion.R")
source("03_MSRun.R")
source("04_DataAnalysis.R")
source("05_Statistics.R")
source("06_Benchmarks.R")
#####################

#####################
## Create list of testing parameters
#####################


# Fasta files:
pathToFasta <- list.files(path = fastaFilePath, full.names = T, pattern = ".fasta")

# Generate a list of parameters to test:
# Note: here we need to define ALL parameters as this will generate the necessary hash
paramGroundTruth <- list(#####   Ground truth generation
  "PathToFasta" = "QC_DataAnalysis/fasta_full_yeast.fasta",
  #"PathToFasta" = pathToFasta,
  "PathToProteinList" = NULL,
  "NumReps" = c(3),
  "NumCond" = 2,
  "FracModProt" = 0,
  "FracModPerProt" = 0,
  "PTMTypes" = NULL,
  "PTMTypesDist" = NULL,
  "PTMTypesMass" = NULL,
  "PTMMultipleLambda" = NULL,
  "ModifiableResidues" = NULL,
  "ModifiableResiduesDistr" = NULL,
  "RemoveNonModFormFrac" = 0
)
paramProteoformAb <- list(
  "QuantNoise" = seq(0.1,0.9,0.5),
  #"QuantNoise" = 0.2,0.7
  "DiffRegFrac" = c(0.1,0.3,0.5),
  "DiffRegMax" = seq(0.5, 2, 0.5),
  "UserInputFoldChanges" = NULL,
  "ThreshNAProteoform" = -100,
  "AbsoluteQuanMean" = 30.5,
  "AbsoluteQuanSD" = 3.6,
  "ThreshNAQuantileProt" = 0.01
)
paramDigest <- list(
  ##### Digestion
  "Enzyme" = "trypsin",
  "PropMissedCleavages" = 0.01, 
  "MaxNumMissedCleavages" = 4,
  "PepMinLength" = 7,
  "PepMaxLength" = 30,
  "LeastAbundantLoss" = 0,
  "EnrichmentLoss" = 0.2,
  "EnrichmentEfficiency" = 1,
  "EnrichmentNonModSignalLoss" = 0,
  "EnrichmentNoise" = 0.2
)
paramMSRun <- list(
  ##### MSRun
  "PercDetectedPep" = seq(0.1,0.5,0.1),
  "PercDetectedVal" = seq(0.1,0.5,0.1),
  "WeightDetectVal" = c(0,0.1,1),
  "MSNoise" = c(0.25, 0.5),
  "WrongIDs" = c(0.01,0.05),
  "WrongLocalizations" = 0.01,
  "MaxNAPerPep" = 1000
)
paramDataAnalysis <- list(
  ##### Data analysis
  "ProtSummarization" = "medpolish",
  "MinUniquePep" = c(1,2),
  "StatPaired" = FALSE
)
#####################

#####################
## Run the sample preparation:
#####################

# all combinations of the parameters
listtogroundtruth <- purrr::cross(compact(paramGroundTruth))

## setting up the entire thing hierarchically.
## Problem: how to parallelize

#Gather always benchmarking data
allBs <- NULL

for (hh in 1:length(listtogroundtruth)) {
  
  ## Running ground thruth generations
  # check whether file with correct parameters exists
  tParam <- listtogroundtruth[[hh]]
  Param <- "none"
  md5 <- digest(tParam, algo = "md5")
  filename <- paste0(resultFilePath,"/outputGroundTruth_",md5,".RData")
  if (file.exists(filename)) {
    load(filename)
  }
  if (all(unlist(Param) != unlist(tParam))) {
    Param <- tParam
    groundTruth <- samplePreparation(parameters = Param)
    save(Param, groundTruth, file = filename)
  }
  gtParam <- Param
  ## quantitative proteoform abundances
  listtoproteoformab <- purrr::cross(compact(paramProteoformAb))
  for (ii in 1:length(listtoproteoformab)) {
    Param <- "none"
    # create combined parameterfile
    tParam <- c(gtParam,listtoproteoformab[[ii]])
    tParam <- c(tParam, list(QuantColnames = paste0("C_",rep(1:tParam$NumCond,each=tParam$NumReps),"_R_", rep(1:tParam$NumReps, tParam$NumCond))))
    # hash code to represent parameter configuration
    md5 <- digest(tParam, algo = "md5")
    filename <- paste0(resultFilePath,"/outputProteoformAb_",md5,".RData")
    if (file.exists(filename)) {
      load(filename)
    }
    # test for exactly same parameters and run if not
    if (!all(unlist(Param) == unlist(tParam))) {
      Param <- tParam
      proteoformAb <- addProteoformAbundance(proteoforms = groundTruth, parameters = Param)
      save(Param, proteoformAb, file = filename)
    }
    pfParam <- Param
    
    ### Digestion
    listtodigestion <- purrr::cross(compact(paramDigest))
    for (jj in 1:length(listtodigestion)) {
      Param <- "none"
      tParam <- c(pfParam,listtodigestion[[jj]])
      md5 <- digest(tParam, algo = "md5")
      filename <- paste0(resultFilePath,"/outputDigestion_",md5,".RData")
      if (file.exists(filename)) {
        load(filename)
      }
      if (!all(unlist(Param) == unlist(tParam))) {
        Param <- tParam
        peptable <- digestGroundTruth(proteoforms = proteoformAb, parameters = c(Param, list(Cores=cores, ClusterType=clusterType)))
        peptable <- digestionProductSummarization(peptides = peptable, parameters = Param)
        BeforeMS <- filterDigestedProt(DigestedProt = peptable, parameters = Param)
        save(Param, BeforeMS, file = filename)
      }
      dgParam <- Param
      
      ### MS run
      listtomsrun <- purrr::cross(compact(paramMSRun))
      for (kk in 1:length(listtomsrun)) {
        Param <- "none"
        tParam <- c(dgParam,listtomsrun[[kk]])
        md5 <- digest(tParam, algo = "md5")
        filename <- paste0(resultFilePath,"/outputMSRun_",md5,".RData")
        if (file.exists(filename)) {
          load(filename)
        }
        if (!all(unlist(Param) == unlist(tParam))) {
          Param <- tParam
          AfterMSRun <- vector(mode = "list")
          for (i in which(sapply(BeforeMS, length) > 0)) {
            AfterMSRun[[length(AfterMSRun) + 1]] <- MSRunSim(Digested = BeforeMS[[i]], parameters = Param)
          }
          names(AfterMSRun) <- names(BeforeMS)[which(sapply(BeforeMS, length) > 0)]
          save(Param, AfterMSRun, file = filename)
        }
        msParam <- Param
        
        ### Protein abundance
        listtodatanalysis <- purrr::cross(compact(paramDataAnalysis))
        for (ll in 1:length(listtodatanalysis)) {
          Param <- "none"
          tParam <- c(msParam,listtodatanalysis[[ll]])
          md5 <- digest(tParam, algo = "md5")
          filename <- paste0(resultFilePath,"/outputDataAnalysis_",md5,".RData")
          if (file.exists(filename)) {
            load(filename)
          }
          if (!all(unlist(Param) == unlist(tParam))) {
            Param <- tParam
            Prots <- proteinSummarisation(peptable = AfterMSRun$NonEnriched, parameters = Param)
            # Don't accept anything below 100 proteins
            if(nrow(Prots) > 99) { 
              # Filter for having at least 1 actual value per protein group and peptide
              Prots <- Prots[rowSums(is.na(Prots[, Param$QuantColnames])) < length(Param$QuantColnames), ]
              allPeps <- as.data.frame(do.call("rbind", AfterMSRun))
              allPeps <- allPeps[rowSums(is.na(allPeps[, Param$QuantColnames])) < length(Param$QuantColnames), ]
              rownames(allPeps) <- paste0("pep", 1:nrow(allPeps))
              Stats <- runPolySTest(Prots, Param, refCond=1, onlyLIMMA=F)
              # much faster with only LIMMA tests   
              StatsPep <- runPolySTest(allPeps, Param, refCond=1, onlyLIMMA=T)
              Benchmarks <- calcBenchmarks(Stats, StatsPep, Param)
              save(Param, Stats, StatsPep, Benchmarks, file = filename)
            } else {
              print("Too few proteins!!!")
              Benchmarks <- Stats <- StatsPep <- NULL
              save(Param, Stats, StatsPep, Benchmarks, file = filename)
            }
          } else if (calcAllBenchmarks) {
            Benchmarks <- calcBenchmarks(Stats, StatsPep, Param)
            save(Param, Stats, StatsPep, Benchmarks, file = filename)
          }
          
          allBs[[md5]] <- list(Benchmarks, Param)
        }
      }
    }
  }
}
cat("###### Finished data set generation \n")

### This part can be used for visualizing and comparing
## Preferably calling an external script

# extracting all benchmarks (sometimes there are more or less per run)
t_allbnames <- NULL
for (i in names(allBs)) {
  t_allbnames <- c(t_allbnames, names(unlist(allBs[[i]][[1]]$globalBMs)))
}
benchNames <- unique(t_allbnames)
BenchMatrix <- data.frame(matrix(NA, ncol=length(benchNames)+length(Param)  , nrow=length(allBs)))
colnames(BenchMatrix) <- c(benchNames, names(Param))
rownames(BenchMatrix) <- names(allBs)

# writing all results and parameters into matrix
for (i in names(allBs)) {
  tglob <- unlist(allBs[[i]][[1]]$globalBMs)
  BenchMatrix[i, names(tglob)] <- tglob
  tpar <- allBs[[i]][[2]]
  BenchMatrix[i, names(tpar)] <- sapply(tpar, function(x) ifelse(length(x)>1, paste0(x,collapse="_"), x))
}

# Visualize roughly
par(mfrow=c(3,3))
# define reference for x-axis
ref <- "WeightDetectVal"
for (obj in benchNames) {
  dat <- BenchMatrix[,obj]
  if (sum(!is.na(dat)) > 0)
    plot(BenchMatrix[,ref], dat, main=obj)
}
par(mfrow=c(1,1))

# ## get all available analyses
# a <- system(paste0("ls ",resultFilePath,"/outputData*"), intern = T)
# BenchMatrix <- data.frame(matrix(NA, ncol=length(benchNames)+length(Param)  , nrow=length(a)))
# colnames(BenchMatrix) <- c(benchNames, names(Param))
# rownames(BenchMatrix) <- names(a)
# 
# 
# for (filename in a) {
#   load(filename)
#   tB <- list(Benchmarks, Param)
#   tglob <- unlist(tB[[1]]$globalBMs)
#   BenchMatrix[filename, names(tglob)] <- tglob
#   tpar <- tB[[2]]
#   BenchMatrix[filename, names(tpar)] <- sapply(tpar, function(x) ifelse(length(x)>1, paste0(x,collapse="_"), x))
#   
# }
# 
# 
# # Visualize roughly
# par(mfrow=c(3,3))
# # define reference for x-axis
# ref <- "DiffRegFrac"
# for (obj in benchNames) {
#   dat <- BenchMatrix[,obj]
#   if (sum(!is.na(dat)) > 0)
#     hist(dat, 100, main=obj)
# }
# par(mfrow=c(1,1))
          