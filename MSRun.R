##################################
## TODO remove afterwards
# Percentage of detected peptides
Param <- list()
Param$PercDetectedPep <- 0.8

# Percentage of detected values (replicate/condition)
Param$PercDetectedVal <- 0.8

# Weights for intensity-dependence of non-detection
Param$WeightDetectVal <- 0

# Wrong identifications
Param$WrongIDs <- 0.01

# Wrong localizations
Param$WrongLocalizations <- 0.01

#############################
Digested <- read.csv("...")

# TODO: what about multiples from different fractions?

# Remove certain percentage of peptides
# Sample number of peptides to be removed
remove <- sample(1:nrow(Digested), Param$PercDetectedPep)
MSRun <- Digested[-remove,]

# Sample number of values to be removed
allVals <- MSRun[,...]

remove <- sample(1:nrow(Digested), Param$PercDetectedPep)
MSRun <- Digested[-remove,]