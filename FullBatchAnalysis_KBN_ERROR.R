################################################################################
#                       QC batch test of the parameters                        #
################################################################################


## only once for installation
 myPaths <- .libPaths()   # get the paths
 #.libPaths("/home/rstudio/R")
 .libPaths("/work/PhosFake Program and results/PhosFake/lib/")
# library("BiocManager")
#BiocManager::install(c("protr","preprocessCore","matrixStats","extraDistr","fdrtool","qvalue","limma","moments"))

 library(purrr)
 library(crayon)
 library(parallel)
 library(digest)
 library(tidyr)
 
#### YOU NEED TO BE IN THE MAIN FOLDER OF THE PHOSFAKE SCRIPTS

####################### Paths and directories
#####################
#Working directory should be PhosFake's main directory.
# File paths
Param <- NULL
fastaFilePath <- "InputData"
resultFilePath <- "KBN_output/KBN_ALBU_BOVIN"
# For parallel computing
cores <- 4
clusterType <- "FORK"
calcAllBenchmarks <- T

try(  dir.create(resultFilePath))


repportName <- paste0(resultFilePath,"/ALBU_Report_BatchAnalysis_", gsub(":|[[:space:]]", "_", format(Sys.time(), "%Y%b%d%X")))

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
  "PathToProteinList" = NULL,
  "NumReps" = c(5), #DEFAULT 3
  "NumCond" = 10,  #DEFAULT 2
  "FracModProt" = 0.2,  #DEFAULT 0.2
  "FracModPerProt" = 0.1,  #DEFAULT 0.1
  "PTMTypes" = "ph", # Phosphorylation
  "PTMTypesDist" = 1, # Frequency distribution of PTM
  "PTMTypesMass" = 72, # Mass difference after the PTM 
  "PTMMultipleLambda" = 0.5, # Lambda value for the poisson distribution used in adjusting probability wights of possible PTM locations
  "ModifiableResidues" = "S", # Which residues are modifiable 
  "ModifiableResiduesDistr" = 1, # Distribution of modifiable residues
  "RemoveNonModFormFrac" = 0 
)
paramProteoformAb <- list(
  "QuantNoise" = c(0.125),
  "DiffRegFrac" = c(0.25, 0.5), 
  "DiffRegMax" = c(3), 
  "UserInputFoldChanges_NumRegProteoforms" = NULL, # rep(2, 2),
  "UserInputFoldChanges_RegulationFC" = NULL , #rep(log2(c(100, 10)), rep(2, 2)), #### use this
  "ThreshNAProteoform" = -100, # remove values below this threshold
  "AbsoluteQuanMean" = 30.5, #DEFAULT
  "AbsoluteQuanSD" = 3.6, #DEFAULT
  "ThreshNAQuantileProt" = 0.0 # Threshold of proteoforms to continue with
)
paramDigest <- list(
  ##### Digestion
  "Enzyme" = "trypsin",
  "PropMissedCleavages" = 0, 
  "MaxNumMissedCleavages" = 0,
  "PepMinLength" = 7, #DEFAULT
  "PepMaxLength" = 30, #DEFAULT
  "LeastAbundantLoss" = 0,
  "EnrichmentLoss" = 0.0,
  "EnrichmentEfficiency" = 0.2,
  "EnrichmentNonModSignalLoss" = 0,
  "EnrichmentNoise" = 0
)
paramMSRun <- list(
  ##### MSRun
  "PercDetectedPep" = 0.1,
  "PercDetectedVal" = 0.9,
  "WeightDetectVal" = 1,
  "MSNoise" = 0.2,
  "WrongIDs" = 0.0,
  "WrongLocalizations" = 0.0,
  "MaxNAPerPep" = 1000
)
paramDataAnalysis <- list(
  ##### Data analysis
  "ProtSummarization" = "medpolish",
  "MinUniquePep" = c(1),
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
              Stats <- runPolySTest(Prots, Param, refCond=1, onlyLIMMA=T)
            # much faster with only LIMMA tests   
              StatsPep <- runPolySTest(allPeps, Param, refCond=1, onlyLIMMA=T)
              Benchmarks <- calcBenchmarks(Stats, StatsPep, Param)
              save(Param, Stats, StatsPep, Benchmarks, file = filename)
            } else {
              print("Too few proteins!!!")
              Benchmarks <- Stats <- StatsPep <- NULL
              save(Param, Stats, StatsPep, Benchmarks, file = filename)
            }
          } else if (calcAllBenchmarks & !is.null(Stats)) {
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

# # Visualize roughly
# par(mfrow=c(3,3))
# # define reference for x-axis
# ref <- "WeightDetectVal"
# for (obj in benchNames) {
#   dat <- BenchMatrix[,obj]
#   if (sum(!is.na(dat)) > 0)
#     plot(BenchMatrix[,ref], dat, main=obj)
# }
# par(mfrow=c(1,1))

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


# need to select from the BeforeMS tibble the name of the sequence and get the control/replicates to plot
##xvaluesBeforeMS <- as.data.frame(BeforeMS[["NonEnriched"]][, c(15:20)])
##xvaluesBeforeMS <- xvaluesBeforeMS - mean(unlist(xvaluesBeforeMS[,1:6]))
##xvaluesBeforeMS <- xvaluesBeforeMS / unlist(xvaluesBeforeMS[which.max(abs(unlist(xvaluesBeforeMS[,1:6])))])
#plot(1:6,xvaluesBeforeMS, xaxt="n" , type = "b", xlab = "Condition # and Replicate #", ylab = "Relative change") 
#matplot(t(xvaluesBeforeMS[1:200,]), type = "l", xaxt="n",, xlab = "Condition # and Replicate #", ylab = "peptide abundance")


condAndRep <- paramGroundTruth[["NumCond"]] * paramGroundTruth[["NumReps"]]+14

# Plot of peptides in the sample before MS run
matplot(t(as.data.frame(BeforeMS[["NonEnriched"]][, c(15:condAndRep)]) - rowMeans(as.data.frame(BeforeMS[["NonEnriched"]][, c(15:20)]), na.rm = T)), 
        type = "b", 
        xaxt="n", 
        xlab = "Condition # and Replicate #", 
        ylab = "Relative abundance", 
        main = paste(length(t(BeforeMS[["NonEnriched"]][,1])), " peptides in samples before MS run"))
axis(1, at=1:6, colnames(BeforeMS[["NonEnriched"]][c(15:condAndRep)]))


##xvaluesAfterMSRun <- as.data.frame(AfterMSRun[["NonEnriched"]][, c(15:20)])
##xvaluesAfterMSRun <- xvaluesAfterMSRun - rowMeans(xvaluesAfterMSRun, na.rm = T)
##xvaluesAfterMSRun <- xvaluesAfterMSRun / unlist(xvaluesAfterMSRun[which.max(abs(unlist(xvaluesAfterMSRun[,1:6])))])
# Plot of peptides identified in MS run
matplot(t(as.data.frame(AfterMSRun[["NonEnriched"]][, c(15:condAndRep)]) - rowMeans(as.data.frame(AfterMSRun[["NonEnriched"]][, c(15:20)]), na.rm = T)),
        type = "b", 
        xaxt="n", 
        xlab = "Condition # and Replicate #", 
        ylab = "Relative abundance", 
        main = paste(length(t(AfterMSRun[["NonEnriched"]][,1])), " peptides identified after MS run"))
axis(1, at=1:6, colnames(BeforeMS[["NonEnriched"]][c(15:condAndRep)]))


AfterMSRunNoNA <- AfterMSRun[["NonEnriched"]] %>% filter(!is.na(Regulation_Amplitude))
AfterMSRunNoNAUnnest <- unnest(AfterMSRunNoNA, c(Peptide, Start, Stop, MC, Accession, Proteoform_ID, Regulation_Amplitude, Regulation_Pattern))

##AfterMSRunNoNA$Regulation_Amplitude <- unlist(AfterMSRunNoNA$Regulation_Amplitude)


#plotting the peptides with ragulated amplitudes
matplot(t(AfterMSRunNoNAUnnest[, c(15:condAndRep)] - rowMeans(AfterMSRunNoNAUnnest[, c(15:20)], na.rm = T)),
        type = "b", 
        xaxt="n", 
        xlab = "Condition # and Replicate #", 
        ylab = "Relative abundance", 
        main = paste(length(t(AfterMSRunNoNAUnnest[,1])), " peptides identified after MS run which are differentially regulated"))
axis(1, at=1:6, colnames(BeforeMS[["NonEnriched"]][c(15:condAndRep)]))

#plotting the peptides with ragulation amplitudes above Cutoff value
cutoffValue <- 2.9
matplot(t(AfterMSRunNoNAUnnest[AfterMSRunNoNAUnnest$Regulation_Amplitude > cutoffValue,][, c(15:condAndRep)] - rowMeans(AfterMSRunNoNAUnnest[AfterMSRunNoNAUnnest$Regulation_Amplitude > cutoffValue,][, c(15:condAndRep)], na.rm = T)),
        type = "b", 
        xaxt="n", 
        xlab = "Condition # and Replicate #", 
        ylab = "Relative abundance", 
        main = paste(length(t(AfterMSRunNoNAUnnest[AfterMSRunNoNAUnnest$Regulation_Amplitude > cutoffValue,][,1])), "peptides identified after MS run with a regulation above", cutoffValue))
axis(1, at=1:6, colnames(BeforeMS[["NonEnriched"]][c(15:condAndRep)]))


NewAfterFiltering <- AfterMSRunNoNAUnnest[rowSums(is.na(AfterMSRunNoNAUnnest[,15:condAndRep])) < 3 , ]
heatmap(as.matrix(NewAfterFiltering[,15:condAndRep]))

NewAfterFilteringAll <- AfterMSRun[[1]][rowSums(is.na(AfterMSRun[[1]][,15:condAndRep])) < 3 , ]
heatmap(as.matrix(NewAfterFilteringAll[,15:condAndRep]))

# make heatmap of diffreg proteins with sidebar which seperates into the peptides which
# actually are differentially regulated

# make single protein plot
# get the accession of the most regulated peptide found and use grep to collect all peptides from the protein accession code
# 
accessionOfInterest <- AfterMSRunNoNAUnnest[which.max(AfterMSRunNoNAUnnest$Regulation_Amplitude),]$Accession
AfterMSRun$NonEnriched[grep(accessionOfInterest, AfterMSRun$NonEnriched$Accession),]
BeforeMS$NonEnriched[grep(accessionOfInterest, BeforeMS$NonEnriched$Accession),]
groundTruth[grep(accessionOfInterest,groundTruth$Accession),]

matplot(t(AfterMSRunNoNA[grep(accessionOfInterest, AfterMSRunNoNA$Accession),][,15:condAndRep] - rowMeans(AfterMSRunNoNA[grep(accessionOfInterest, AfterMSRunNoNA$Accession),][,15:condAndRep], na.rm = T)),
        type = "b", 
        xaxt="n", 
        xlab = "Condition # and Replicate #", 
        ylab = "Relative abundance", 
        main = paste("The most differentially regulated protein is", accessionOfInterest, "\n which has", length(AfterMSRunNoNA[grep(accessionOfInterest, AfterMSRunNoNA$Accession),][,1]), "peptides"))
axis(1, at=1:6, colnames(BeforeMS[["NonEnriched"]][c(15:condAndRep)]))
