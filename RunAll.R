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
#####################
