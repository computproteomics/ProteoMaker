################################################################################
#                     TO TUN THE ENTIRE PHOSFAKE PIPELINE                      #
################################################################################

#####################
## Load parameters
#####################
source("Parameter.R")
#####################

#####################
## Run the sample preparation simulation:
#####################
source("01_GenerateGroundTruth.R")
# Create the initial list of proteoforms:
proteoforms <- samplePreparation(fasta.path = Param$PathToFasta, parameters = Param)
# Create the full structure of proteoforms along with abundances and ground truth expression patterns:
GroundTruth <- addProteoformAbundance(proteoforms = proteoforms, parameters = Param)
rm(proteoforms)
# Save GroundTruth for analysis:
save(GroundTruth, file = "RData/GroundTruthAbs.RData")
#####################

#####################
## Digestion and sample peptide enrichment:
#####################
source("02_Digestion.R")
# Digest all the proteoforms and get peptide table:
peptable <- DigestGroundTruth(GroundTruth = GroundTruth, parameters = Param)
# Save peptable before filter for analysis:
save(peptable, file = "RData/AllPep.RData")
peptable <- mapQuanToDigestionProd(DigestedProt = peptable)
BeforeMS <- filterDigestedProt(peptable, Param)
# Save peptable before in silico MS run:
save(BeforeMS, file = "RData/BeforeMS.RData")
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
save(AfterMSRun, file = "RData/AfterMSRun.RData")
#####################