################################################################################
#                     TO TUN THE ENTIRE PHOSFAKE PIPELINE                      #
################################################################################

#####################
## Load parameters
#####################
source("Parameter.R")
#####################

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

Stats <- runPolySTest(Prots, Param, refCond=1)

# Optional as it takes quite some time
# allPeps <- do.call("rbind", AfterMSRun)
# StatsPep <- runPolySTest(allPeps, Param, refCond=1)

Stats <- Stats[rowSums(is.na(Stats[, Param$QuantColnames])) < length(Param$QuantColnames), ]

calcROC <- function (columnName, groundtruthColumn="min1Reg") {
FPs <- TPs <- 0
FNs <- sum(Stats$min1Reg)
TNs <- nrow(Stats) - FNs
ROC <- NULL
FDR <- NULL
AUC <- 0
oldFDR <- -1
oldFPR <- 0
for (i in order(Stats[,columnName])) {
  if(Stats[i, groundtruthColumn]) {
        TPs <- TPs + 1
        FNs <- FNs -1
  } else {
    FPs <- FPs + 1
    TNs <- TNs - 1
  }
  currFDR <- Stats[i, columnName]  
  if (is.na(currFDR)) currFDR <- 1
  if (currFDR != oldFDR) {
    ROC <- rbind(ROC, c(FPs/(FPs+TNs), TPs/(TPs+FNs)))
    FDR <- rbind(FDR, c(currFDR, FPs/(FPs + TPs)))
    AUC <- AUC + (FPs/(FPs+TNs) - oldFPR) *  TPs/(TPs+FNs)
  } else {
    ROC[nrow(ROC), ] <- c(FPs/(FPs+TNs), TPs/(TPs+FNs))
    FDR[nrow(FDR), ] <- c(currFDR, FPs/(FPs + TPs))    
    AUC <- AUC + (FPs/(FPs+TNs) - oldFPR) *  TPs/(TPs+FNs)
  }

  oldFDR <- currFDR
  oldFPR <- FPs/(FPs+TNs)
}
cbind(ROC, FDR, AUC)
}

ROC <-  calcROC("FDR PolySTest 2 vs 1")
plot(ROC[,1], ROC[,2], type="s")
lines(ROC[,3], ROC[,4], type="s")
ROC2 <-  calcROC("FDR t-test 2 vs 1")
lines(ROC2[,1], ROC2[,2], col=2, type="s")
lines(ROC2[,3], ROC2[,4], type="s", col=2)
ROC3 <-  calcROC("FDR limma 2 vs 1")
lines(ROC3[,1], ROC3[,2], col=3, type="s")
lines(ROC3[,3], ROC3[,4], type="s", col=3)
ROC4 <-  calcROC("FDR permutation test 2 vs 1")
lines(ROC4[,1], ROC4[,2], col=4, type="s")
lines(ROC4[,3], ROC4[,4], type="s", col=4)
ROC5 <-  calcROC("FDR rank products 2 vs 1")
lines(ROC5[,1], ROC5[,2], col=5, type="s")
lines(ROC5[,3], ROC5[,4], type="s", col=5)
ROC6 <-  calcROC("FDR Miss test 2 vs 1")
lines(ROC6[,1], ROC6[,2], col=6, type="s")
lines(ROC6[,3], ROC6[,4], type="s", col=6)
legend("bottomright", lwd=rep(1,6), col=1:6, cex=0.7, legend = c(paste("PolySTest", ROC[1,5]), paste("t-test", ROC2[1,5]), paste("LIMMA", ROC3[1,5]),
                                           paste("Permutation", ROC4[1,5]),paste("Rank Products", ROC5[1,5]),paste("Miss test", ROC6[1,5])))
  