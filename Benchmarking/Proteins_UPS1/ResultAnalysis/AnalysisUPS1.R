################################################################################
#                     Analysis of PhosFake mimicking UPS1                      #
################################################################################

#####################
## Load parameters
#####################
source("../Parameters/Parameter00.R")
name <- "UPS1_00"
pathToRes <- "../RData/"
#####################


#####################
## Load data
#####################
load(file = paste0(pathToRes, name, "/GroundTruth.RData"))
load(file = paste0(pathToRes, name, "/BeforeMS.RData"))
peptable <- BeforeMS$NonEnriched
#####################


#####################
## Ground truth data set:
#####################
cat("Number of protein accessions in the ground truth data set:\n")
print(length(unique(GroundTruth$Accession)))

cat("Number of proteoforms in the ground truth data set:\n")
print(nrow(GroundTruth))

cat("Number of missingvalues:\n")
print(sum(is.na(GroundTruth[,grepl("^C_", names(GroundTruth))])))
#####################


#####################
## Statistical analysis of the proteoforms. 
#####################

#####################


sessionInfo()