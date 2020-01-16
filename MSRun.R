library(stringr)
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


# residues for PTM type and relative distribution to set them
Param$ModifiableResidues <- list()
Param$ModifiableResiduesDistr <- list()


Param$PTMTypes <- c("ph","ox")

Param$ModifiableResidues[["ph"]] <- c("S","T","Y")
Param$ModifiableResidues[["ox"]] <- c("M")
Param$ModifiableResiduesDistr[["ph"]] <- c(0.86,0.13,0.01)
Param$ModifiableResiduesDistr[["ox"]] <- c(1)



#############################

### will be removed ###
load("data/expDataFrame.RData")
Digested <- finalDF
Digested[, quant_colnames] <- log2 (Digested[, quant_colnames])
######################
# TODO: what about multiples from different fractions?

# Remove certain percentage of peptides
# Sample number of peptides to be removed
remove <- sample(1:nrow(Digested), size = (1-Param$PercDetectedPep)*nrow(Digested))
MSRun <- Digested[-remove,]

# Sample number of values to be remo  ved
allVals <- as.vector(unlist(MSRun[,quant_colnames]))
par(mfrow=c(1,2))
hist(allVals,100)
remove <- sample(1:length(allVals), size=(1-Param$PercDetectedVal)*length(allVals), prob = (rank(-allVals)/length(allVals))  ^  Param$WeightDetectVal)
allVals[remove] <- NA
hist(allVals,100)
par(mfrow=c(1,1))
MSRun[,quant_colnames] <- allVals

# shuffle identifications
shuffle <- sample(1:nrow(MSRun), Param$WrongIDs*nrow(MSRun))
for (i in 1:floor(length(shuffle)/2)) {
  temp <-  MSRun[shuffle[(i-1)*2+1], quant_colnames]
   MSRun[shuffle[(i-1)*2+1], quant_colnames] <- MSRun[shuffle[i*2],quant_colnames]
   MSRun[shuffle[i*2],quant_colnames] <- temp
}
if (length(shuffle)%%2 == 1) {
  temp <-  MSRun[shuffle[1], quant_colnames]
  MSRun[shuffle[1], quant_colnames] <- MSRun[shuffle[length(shuffle)],quant_colnames]
  MSRun[shuffle[length(shuffle)],quant_colnames] <- temp
}

MSRun$WrongID <- F
MSRun$WrongID[shuffle] <- T

## false localizations

# loop over PTM types
for (mod in Param$PTMTypes) {
  isModified <- which(sapply(MSRun$PTMType, function(x) any(x == mod,na.rm=T)))
  modified <- MSRun[isModified, c("PepSequence", "PTMPos", "PTMType")]
  # count number of modifiable residues
    tmpModPosCount <- str_locate_all(modified$PepSequence, paste0(Param$ModifiableResidues[[mod]],collapse="|"))
    ModifiableCount <- sapply(tmpModPosCount, function(x) length(x[,1]))
    ModCount <- sapply(modified$PTMType, function(x) sum(x == mod))
    # positions in original table where modified peptide can get mislocated PTM
    CanMisLoc <- sum(ModifiableCount > ModCount)
    # number of modified peptides where we change localization
    NumForMisLoc <- round(Param$WrongLocalizations*CanMisLoc)
    print(paste0(NumForMisLoc))
    
    # mislocate NOW
    for (ind in sample(isModified[ModifiableCount > ModCount], NumForMisLoc)) {
      curr_pep <- MSRun[ind,]
      str_locate_all(curr_pep$Sequence, paste0(Param$ModifiableResidues[[mod]],collapse="|"))[[1]][,1]
      PTMpos <- curr_pep$PTMPos - curr_pep$PepStart
    }
}


