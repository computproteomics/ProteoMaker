################################################################################
#                     TO TUN THE ENTIRE PHOSFAKE PIPELINE                      #
  ################################################################################
  
  #####################
  ## Load parameters
  #####################
  source("Parameter.R")
  #####################
  

## TEMPORARY
allBs <- NULL
for (p in seq(0,10,0.2)) {
  Param$FracModPerProt <- i
####


  #####################
  ## Option to save data:
  #####################
  if (!exists("saveData")) {
    saveData <- FALSE
  }
  #####################
  
  #####################
  ## Run the sample preparation simulation:
  #####################
  source("01_GenerateGroundTruth.R")
  # Create the initial list of proteoforms:
  proteoforms <- samplePreparation(parameters = Param)
  # Create the full structure of proteoforms along with abundances and ground truth expression patterns:
  GroundTruth <- addProteoformAbundance(proteoforms = proteoforms, parameters = Param)
  rm(proteoforms)
  # Save GroundTruth for analysis:
  if (saveData) {
    save(GroundTruth, file = "RData/GroundTruthAbs.RData")
  }
#####################

#####################
## Digestion and sample peptide enrichment:
#####################
  source("02_Digestion.R")
# Digest all the proteoforms and get peptide table:
peptable <- digestGroundTruth(proteoforms = GroundTruth, parameters = Param)
  # Save peptable before filter for analysis:
if (saveData) {
  save(peptable, file = "RData/AllPep.RData")
}
peptable <- digestionProductSummarization(peptides = peptable, parameters = Param)
BeforeMS <- filterDigestedProt(DigestedProt = peptable, parameters = Param)
# Save peptable before in silico MS run:
if (saveData) {
  save(BeforeMS, file = "RData/BeforeMS.RData")
}
#####################


#####################
## Simulate MS analysis
#####################
source("03_MSRun.R")
AfterMSRun <- vector(mode = "list")
for (i in which(sapply(BeforeMS, length) > 0)) {
  AfterMSRun[[length(AfterMSRun) + 1]] <- MSRunSim(Digested = BeforeMS[[i]], parameters = Param)
}
names(AfterMSRun) <- names(BeforeMS)[which(sapply(BeforeMS, length) > 0)]
# Save final peptable tables for analysis:
if (saveData) {
  save(AfterMSRun, file = "RData/AfterMSRun.RData")
}
#####################


# #####################
# ## Protein abundance calculation
# #####################
source("04_DataAnalysis.R")
Prots <- proteinSummarisation(peptable = AfterMSRun$NonEnriched, parameters = Param)

# #####################
# ## Statistical testing
# ##################### 
source("05_Statistics.R")
  
Stats <- runPolySTest(Prots, Param, refCond=1, onlyLIMMA=F)

allPeps <- as.data.frame(do.call("rbind", AfterMSRun))
# could be dangerous when the same peptide appears in different fractions 
rownames(allPeps) <- paste0("pep", 1:nrow(allPeps))
# much faster with only LIMMA tests   
StatsPep <- runPolySTest(allPeps, Param, refCond=1, onlyLIMMA=T)


# #####################
# ## Collecting benchmarks    
# #####################   
source("06_Benchmarks.R")

# Filter for having at least 1 actual value per protein group and peptide
Stats <- Stats[rowSums(is.na(Stats[, Param$QuantColnames])) < length(Param$QuantColnames), ]
StatsPep <- StatsPep[rowSums(is.na(StatsPep[, Param$QuantColnames])) < length(Param$QuantColnames), ]

Benchmarks <- calcBenchmarks(Stats, StatsPep, Param)
    allBs[[p]] <- Benchmarks
}