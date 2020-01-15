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


# Number of conditions
Param$NumCond <- 3
# Number of replicates
Param$NumReps <- 4

quant_colnames <- paste0("C_",rep(1:Param$NumCond,each=Param$NumReps),"_R_", rep(1:Param$NumReps, Param$NumCond))

#############################

### will be removed ###
load("data/expDataFrame.RData")
Digested <- formattedDFcast[,c(1:8,11:20)]
colnames(Digested)[5:16] <- paste0("C_",rep(1:3,each=4),"_R_", rep(1:4, 3))
######################
# TODO: what about multiples from different fractions?

# Remove certain percentage of peptides
# Sample number of peptides to be removed
remove <- sample(1:nrow(Digested), size = Param$PercDetectedPep*nrow(Digested))
MSRun <- Digested[-remove,]

# Sample number of values to be remo  ved
allVals <- as.vector(unlist(MSRun[,quant_colnames]))

hist(allVals)
remove <- sample(1:length(allVals), size=Param$PercDetectedVal*length(allVals), prob = rank(allVals)/length(allVals)  ^  Param$WeightDetectVal)
allVals[remove] <- NA
hist(allVals)

remove <- sample(1:nrow(Digested), Param$PercDetectedPep)
MSRun <- Digested[-remove,]