################################################################################
#                     TO TUN THE ENTIRE PHOSFAKE PIPELINE                      #
################################################################################

#####################
## Load parameters
#####################
source("../Parameters/Parameter04.R")
pathToRes <- "../RData/UPS1_04/"
#####################

#####################
## Run the sample preparation simulation:
#####################
source("../../../01_GenerateGroundTruth.R")
# Create the initial list of proteoforms:
proteoforms <- samplePreparation(fasta.path = Param$PathToFasta, parameters = Param)
# Create the full structure of proteoforms along with abundances and ground truth expression patterns:
GroundTruth <- addProteoformAbundance(proteoforms = proteoforms, parameters = Param)
rm(proteoforms)
# Save GroundTruth for analysis:
save(GroundTruth, file = paste0(pathToRes, "GroundTruth.RData"))
#####################

#####################
## Digestion and sample peptide enrichment:
#####################
source("../../../02_Digestion.R")
# Digest all the proteoforms and get peptide table:
peptable <- DigestGroundTruth(GroundTruth = GroundTruth, parameters = Param)
save(peptable, file = paste0(pathToRes, "AllPep.RData"))
peptable <- mapQuanToDigestionProd(DigestedProt = peptable)
BeforeMS <- filterDigestedProt(peptable, Param)
# # Save peptable before filter for analysis:
save(BeforeMS, file = paste0(pathToRes, "BeforeMS.RData"))
#####################

sessionInfo()
