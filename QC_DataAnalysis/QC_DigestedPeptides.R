################################################################################
#                     QC of the Ground truth proteoforms                       #
################################################################################

#####################
## Load parameters
#####################
setwd("/Users/rpk349/Documents/Boulot/GitRepo/PhosFake")
source("Parameter.R")
#####################


#####################
## Load ground truth data (needs to be saved after beeing generated in RunAll.R)
#####################
load("RData/peptides.RData")
#####################

#####################
## Peptides 
#####################
peptable <- tosave
boxplot(peptable[,grepl("^C_", names(peptable))])
#####################
