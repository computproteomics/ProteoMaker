################################################################################
#                               IN SILICO MS RUN                               #
################################################################################

library(stringr)

MSRunSim <- function(Digested, parameters) {
  
  # TODO: what about multiples from different fractions?
  
  # Add noise to MS analysis
  # random noise
  matnoise <- rnorm(n = nrow(Digested), mean = 0, sd = parameters$MSNoise)
  for (i in seq_along(parameters$QuantColnames)) {
    
    name <- parameters$QuantColnames[i]
    Digested[,name] = Digested[,name] + matnoise[i]
    
  }
  
  # Remove certain percentage of peptides
  # Sample number of peptides to be removed
  remove <- sample(1:nrow(Digested), size = (1-parameters$PercDetectedPep)*nrow(Digested))
  MSRun <- Digested[-remove,]
  
  # Sample number of values to be removed
  allVals <- as.vector(unlist(MSRun[,parameters$QuantColnames]))
  # par(mfrow=c(1,2))
  # hist(allVals,100)
  if(parameters$WeightDetectVal > 0) {
    myprob <-  (rank(-allVals)/length(allVals))  ^  parameters$WeightDetectVal
    # TODO: sample is too slow with "replace = F". So I iterate to use replace = T. 
    valtoiter <- 1:length(allVals)
    remove <- sample(valtoiter, size=(1-parameters$PercDetectedVal)*length(valtoiter), prob = myprob, replace = T)
    while (sum(duplicated(remove)) > 0 | is.null(remove)) {
      valtoiter <- setdiff(valtoiter, unique(remove))
      myprobiter <- myprob[valtoiter]
      remove[duplicated(remove)] <- sample(valtoiter, size=sum(duplicated(remove)), prob = myprobiter, replace = T)
    }
  } else {
    cat("WARNING: missing values at peptide level (due to in silico MS detection) are determined at random.\n")
    cat("         Change the parameter \'WeightDetectVal\' to a value > 0 to add weight to lowest intensities.")
    remove <- sample(1:length(allVals), size=(1-parameters$PercDetectedVal)*length(allVals))
  }

  
  allVals[remove] <- NA
  # hist(allVals,100)
  # par(mfrow=c(1,1))
  MSRun[,parameters$QuantColnames] <- allVals
  
  # shuffle identifications
  shuffle <- sample(1:nrow(MSRun), parameters$WrongIDs*nrow(MSRun))
  for (i in 1:floor(length(shuffle)/2)) {
    temp <-  MSRun[shuffle[(i-1)*2+1], parameters$QuantColnames]
    MSRun[shuffle[(i-1)*2+1], parameters$QuantColnames] <- MSRun[shuffle[i*2],parameters$QuantColnames]
    MSRun[shuffle[i*2],parameters$QuantColnames] <- temp
  }
  if (length(shuffle)%%2 == 1) {
    temp <-  MSRun[shuffle[1], parameters$QuantColnames]
    MSRun[shuffle[1], parameters$QuantColnames] <- MSRun[shuffle[length(shuffle)],parameters$QuantColnames]
    MSRun[shuffle[length(shuffle)],parameters$QuantColnames] <- temp
  }
  
  MSRun$WrongID <- F
  MSRun$WrongID[shuffle] <- T
  
  ## false localizations
  # loop over PTM types
  MSRun$IsMisLocated <- F
  if (parameters$FracModProt > 0) {
    for (mod in 1:length(parameters$PTMTypes)) {
      isModified <- which(sapply(MSRun$PTMType, function(x) any(x == parameters$PTMTypes[mod],na.rm=T)))
      modified <- MSRun[isModified, c("Sequence", "PTMPos", "PTMType")]
      # count number of modifiable residues
      tmpModPosCount <- str_locate_all(modified$Sequence, paste0(parameters$ModifiableResidues[[mod]],collapse="|"))
      ModifiableCount <- sapply(tmpModPosCount, function(x) length(x[,1]))
      ModCount <- sapply(modified$PTMType, function(x) sum(x == parameters$PTMTypes[mod]))
      # positions in original table where modified peptide can get mislocated PTM
      CanMisLoc <- sum(ModifiableCount > ModCount)
      # number of modified peptides where we change localization
      NumForMisLoc <- round(parameters$WrongLocalizations*CanMisLoc)
      print(paste0(NumForMisLoc))
      
      # mislocate NOW
      if (NumForMisLoc > 0) {
        for (ind in sample(isModified[ModifiableCount > ModCount], NumForMisLoc)) {
          curr_pep <- MSRun[ind,]
          available_pos <- str_locate_all(curr_pep$Sequence, paste0(parameters$ModifiableResidues[[mod]],collapse="|"))[[1]][,1]
          PTMpos <- curr_pep$PTMPos[[1]][curr_pep$PTMType[[1]] == parameters$PTMTypes[mod]]
          newPos <- sample(available_pos[!(available_pos %in% PTMpos)],1)
          tPTMvec <- curr_pep$PTMPos[[1]]
          if(length(PTMpos) != 1){PTMpos <- sample(PTMpos,1)}
          tPTMvec[tPTMvec == PTMpos] <- newPos
          curr_pep$PTMPos <- list(tPTMvec)
          curr_pep$IsMisLocated <- T
          MSRun[ind,] <- curr_pep
        }
      }
    }
  }
  
  return(MSRun)
}