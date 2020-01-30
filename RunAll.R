################################################################################
#                     TO TUN THE ENTIRE PHOSFAKE PIPELINE                      #
################################################################################

## Load parameters
#####################
source("Parameter.R")
#####################


## Run the sample preparation simulation:
#####################
source("SamplePrep.R")
#####################

## Sample preparation
#####################
#Use just to create proteoforms (modified and unmodified)
# proteoforms <- samplePreparation(fasta.path = "fasta.fasta", parameters = Param)
#Use to create the full structure of proteoforms along with abundances and ground truth expression patterns.

# #For example file (50 proteins)
# GroundTruth = addProteoformAbundance(proteoforms = samplePreparation(fasta.path = "fasta_example.fasta", parameters = Param), parameters = Param)
#For whole human proteome
GroundTruth = addProteoformAbundance(proteoforms = samplePreparation(fasta.path = "fasta_full_human.fasta", parameters = Param), parameters = Param)

#####################
