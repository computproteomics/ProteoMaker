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

