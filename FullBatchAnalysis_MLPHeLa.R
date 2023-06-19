################################################################################
#                       QC batch test of the parameters                        #
################################################################################

library(purrr)
library(crayon)
library(parallel)
library(digest)

# # only once for installation
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
  "PathToFasta" = "QC_DataAnalysis/fasta_full_human.fasta",
  "PathToProteinList" = NULL,
  "NumReps" = 1,
  "NumCond" = 1,
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
  "QuantNoise" = seq(0.01,0.6,0.1),
  "DiffRegFrac" = 0,
  "DiffRegMax" = 0,
  "UserInputFoldChanges" = NULL,
  "ThreshNAProteoform" = -100,
  "AbsoluteQuanMean" = 30.5, # default value: 30.5
  "AbsoluteQuanSD" = seq(2.5,3.5,0.5), # default value: 3.6
  "ThreshNAQuantileProt" = seq(0.01,0.6, 0.1)
)
paramDigest <- list(
  ##### Digestion
  "Enzyme" = "trypsin",
  "PropMissedCleavages" = 0, 
  "MaxNumMissedCleavages" = 0,
  "PepMinLength" = 7,
  "PepMaxLength" = 30,
  "LeastAbundantLoss" = seq(0.03, 0.6, 0.1),
  "EnrichmentLoss" = 0,
  "EnrichmentEfficiency" = 1,
  "EnrichmentNonModSignalLoss" = 0,
  "EnrichmentNoise" = 0
)
paramMSRun <- list(
  ##### MSRun
  "PercDetectedPep" = seq(0.4, 0.7, 0.1),
  "PercDetectedVal" = seq(0.4, 0.7, 0.1),
  "WeightDetectVal" = c(0,1,2,3,4,5,10,20,50,100),
  "MSNoise" = seq(0.01, 0.6, 0.1),
  "WrongIDs" = 0.01,
  "WrongLocalizations" = 0.01,
  "MaxNAPerPep" = 1000
)
paramDataAnalysis <- list(
  ##### Data analysis
  "ProtSummarization" = "medpolish",
  "MinUniquePep" = 1,
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
        # msParam <- Param
        # 
        # ### Protein abundance
        # listtodatanalysis <- purrr::cross(compact(paramDataAnalysis))
        # for (ll in 1:length(listtodatanalysis)) {
        #   Param <- "none"
        #   tParam <- c(msParam,listtodatanalysis[[ll]])
        #   md5 <- digest(tParam, algo = "md5")
        #   filename <- paste0(resultFilePath,"/outputDataAnalysis_",md5,".RData")
        #   if (file.exists(filename)) {
        #     load(filename)
        #   }
        #   if (!all(unlist(Param) == unlist(tParam))) {
        #     Param <- tParam
        #     Prots <- proteinSummarisation(peptable = AfterMSRun$NonEnriched, parameters = Param)
        #     # Don't accept anything below 100 proteins
        #     if(nrow(Prots) > 99) {
        #       # Filter for having at least 1 actual value per protein group and peptide
        #       Prots <- Prots[rowSums(is.na(Prots[, Param$QuantColnames])) < length(Param$QuantColnames), ]
        #       allPeps <- as.data.frame(do.call("rbind", AfterMSRun))
        #       allPeps <- allPeps[rowSums(is.na(allPeps[, Param$QuantColnames])) < length(Param$QuantColnames), ]
        #       rownames(allPeps) <- paste0("pep", 1:nrow(allPeps))
        #       Stats <- runPolySTest(Prots, Param, refCond=1, onlyLIMMA=T)
        #       # much faster with only LIMMA tests
        #       StatsPep <- runPolySTest(allPeps, Param, refCond=1, onlyLIMMA=T)
        #       Benchmarks <- calcBenchmarks(Stats, StatsPep, Param)
        #       save(Param, Stats, StatsPep, Benchmarks, file = filename)
        #     } else {
        #       print("Too few proteins!!!")
        #       Benchmarks <- Stats <- StatsPep <- NULL
        #       save(Param, Stats, StatsPep, Benchmarks, file = filename)
        #     }
        #   } else if (calcAllBenchmarks & !is.null(Stats)) {
        #     Benchmarks <- calcBenchmarks(Stats, StatsPep, Param)
        #     save(Param, Stats, StatsPep, Benchmarks, file = filename)
        #   }
        #   
        #   allBs[[md5]] <- list(Param)
        # }
        allBs[[md5]] <- list(Param)
      }
    }
  }
}

cat("###### Finished data set generation \n")

### This part can be used for visualizing and comparing
## Preferably calling an external script

# extracting all benchmarks (sometimes there are more or less per run)
BenchMatrix <- data.table::rbindlist(unlist(allBs, recursive = F))
BenchMatrix <- as.data.frame(BenchMatrix)
BenchMatrix$FileName <- names(allBs)

# >> MLP: I save the table for future analysis in comparison to HeLa runs:

repportName <- paste0("QC_DataAnalysis/HeLaLikeSimulation/Parameters_", gsub(":|[[:space:]]", "_", format(Sys.time(), "%Y%b%d%X")), ".txt")

write.table(BenchMatrix, file = repportName, sep = "\t", row.names = F, quote = F)
cat("-> Parameter saved:", repportName, "\n")

# # >> MLP: Open old repports:
# 
# BenchMatrix <- read.table(file = "QC_DataAnalysis/HeLaLikeSimulation/Parameters_2020Jun0208_10_57_AM.txt",
#                           sep = "\t",
#                           header = T)

## Benchmarks:

# wrapper for calculating metrics related to peptides only, with no stat:
calcBenchmarksNoStat <- function(AfterMSRun)  {
  
  Benchmarks <- NULL
  
  ## Remove peptides with only missing values:
  AfterMSRun <- AfterMSRun[!is.na(AfterMSRun$C_1_R_1),]
  
  # global means 1 number per metric
  globalBMs <- list(
    # Peptide level
    numPeptides=0, numProteins=0, propUniquePep=0, propSharedPep=0, percMissingPep=0,
    propMisCleavedPeps=0,skewnessPeps=0, kurtosisPeps=0, sdPeps=0, numModPeptides=0)  
  
  #### Calculating peptide numbers
  globalBMs["numPeptides"] <- nrow(AfterMSRun)
  globalBMs["numProteins"] <- length(unique(unlist(AfterMSRun$Accession)))
  globalBMs["propUniquePep"] <-  sum(sapply(AfterMSRun$Accession, function(x) length(x) == 1)) / nrow(AfterMSRun)
  globalBMs["propSharedPep"] <-  sum(sapply(AfterMSRun$Accession, function(x) length(x) > 1)) / nrow(AfterMSRun)
  globalBMs["percMissingPep"] <- sum(is.na(unlist(AfterMSRun[,grepl("C.+_R.+", names(AfterMSRun))]))) / length(unlist(AfterMSRun[,grepl("C.+_R.+", names(AfterMSRun))])) * 100
  
  # miscleavages
  globalBMs["propMisCleavedPeps"] <- sum(sapply(AfterMSRun$MC, unique)) / nrow(AfterMSRun)
  
  # checking properties of distribution
  globalBMs$skewnessPeps <- skewness(unlist(AfterMSRun[,grepl("C.+_R.+", names(AfterMSRun))]), na.rm=T)
  globalBMs$kurtosisPeps <- kurtosis(unlist(AfterMSRun[,grepl("C.+_R.+", names(AfterMSRun))]), na.rm=T)
  globalBMs$sdPeps <- sd(unlist(AfterMSRun[,grepl("C.+_R.+", names(AfterMSRun))]), na.rm=T)
  
  globalBMs$numModPeptides <- sum(AfterMSRun$PTMType != "NULL")
  
  Benchmarks$globalBMs <- globalBMs
  return(Benchmarks)
  
}

nr <- nrow(BenchMatrix)
matres <- matrix(nrow = nr, ncol = 10)
# startat <- 1
startat <- 7108
matres <- matrix(nrow = (nr - startat + 1), ncol = 10)

for (i in startat:nr) { # testing
  cat("\n#### Calculating benchmarks", i, "( of", nr, ") ###\n")
  filename <- paste0(resultFilePath,"/outputMSRun_",BenchMatrix$FileName[i],".RData")
  load(filename)
  matres[i,] <- unlist(calcBenchmarksNoStat(AfterMSRun = AfterMSRun[[1]])[[1]])
}
colnames(matres) <- names(unlist(calcBenchmarksNoStat(AfterMSRun = AfterMSRun[[1]])[[1]]))

## Add Benchmarks to the parameter table:

tab_benchmarksHeLa <- as.data.frame(matres)
tab_benchmarksHeLa$FileName <- BenchMatrix$FileName[startat:nr]

# >> MLP: I save the table for output analysis:

repportName <- paste0("QC_DataAnalysis/HeLaLikeSimulation/Results_", gsub(":|[[:space:]]", "_", format(Sys.time(), "%Y%b%d%X")), ".txt")

write.table(tab_benchmarksHeLa, file = repportName, sep = "\t", row.names = F, quote = F)
cat("-> Results saved:", repportName, "\n")


#>>> REMOVE FILES:

toremove <- bm$FileName[bm$AbsoluteQuanSD < 2.4 | bm$AbsoluteQuanSD > 3.6]
myfiles <- list.files(resultFilePath, full.names = T, recursive = T)

for (na in toremove) { 
  filename <- myfiles[grepl(na,myfiles)]
  sapply(filename, file.remove)
}
