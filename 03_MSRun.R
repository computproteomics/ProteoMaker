################################################################################
#                               IN SILICO MS RUN                               #
################################################################################

library(stringr)
library(crayon)

#####################
## Function that simulates MS run analysis.
## - Adds random noise to each sample due to MS instrument.
## - Removes a specific proportion of peptides due to detection limitations.
## - Removes a percentage of non-missing intensities based on probability weights depending on intensity values.
## - Introduces peptide identification false discovery rate.
## - Adds PTMs false localization rate for multiple PTM types.
## - Filters out peptides based on a maximum number of missing values threshold.
#####################
MSRunSim <- function(Digested, parameters) {
  
  # TODO: what about multiples from different fractions?
  
  cat("#MS RUN SIMULATION - Start\n\n")
  cat(" + Noise addition:\n")
  cat("  - The MS noise standard deviation is", parameters$MSNoise, ".\n")

  # Introducing random noise to MS analysis due to MS instrument.
  matnoise <- rnorm(n = nrow(Digested) * length(parameters$QuantColnames), mean = 0, sd = parameters$MSNoise)

#  for (i in seq_along(parameters$QuantColnames)) {

#    Digested[ ,parameters$QuantColnames[i]] <- Digested[ ,parameters$QuantColnames[i]] + matnoise[i]
    Digested[ ,parameters$QuantColnames] <- Digested[ ,parameters$QuantColnames] + matnoise[i]

#  }
  
  cat("  - Noise added to all samples!\n\n")
  cat(" + Detection limits:\n")
  cat("  - The percentage of remaining peptides is", parameters$PercDetectedPep*100, "%.\n")
  
  # Sample a percentage of random peptides to be removed.
  if(parameters$PercDetectedPep != 1) {
    
    remove <- sample(1:nrow(Digested), size = (1-parameters$PercDetectedPep)*nrow(Digested))
    MSRun <- Digested[-remove, ]
    cat("  - A total of", (1-parameters$PercDetectedPep)*nrow(Digested), "peptides is removed.\n\n")
  
  } else {
    
    MSRun <- Digested
    cat("  - No peptides were removed.\n\n")
    
  }
  
  # Sample a random percentage of intensities to be removed.
  allVals <- as.vector(unlist(MSRun[ ,parameters$QuantColnames]))
  cat(" + Removing peptide intensities:\n")
  cat("  - The percentage of remaining intensities is", parameters$PercDetectedVal*100, "% and the probabilities of selection depends on intensity values by", parameters$WeightDetectVal, ".\n")
  myprob <-  (rank(-allVals)/length(allVals)) ^ parameters$WeightDetectVal

  #The probability for existing NA values from the above will be 1.
  #Thus, replacing probabilities corresponding to existing NAs to 0.
  myprob[is.na(allVals)] <- 0
  
  #This is the fastest way to sample with replacement. 
  #Reference: https://doi.org/10.1016/j.ipl.2005.11.003
  remove <- order(runif(length(allVals)) ^ (1/myprob), decreasing = T)[1:((1-parameters$PercDetectedVal)*length(allVals))]
  
  if(parameters$WeightDetectVal == 0){
    
    cat(crayon::red("WARNING: missing values at peptide level (due to in silico MS detection) are determined at random.\n"))
    cat(crayon::red("         Change the parameter \'WeightDetectVal\' to a value > 0 to add weight to lowest intensities."))
    
  }
  
  allVals[remove] <- NA
  MSRun[ ,parameters$QuantColnames] <- allVals
  cat("  - A total of", length(remove), "intensities is removed.\n\n")
  
  # Shuffle the intensities of randomly selected peptides, to express the wrong identification.
  shuffle <- order(runif(nrow(MSRun)), decreasing = T)[1:(parameters$WrongIDs*nrow(MSRun))]
  cat(" + Addition of false identification:\n")
  cat("  - FDR selected is", parameters$WrongIDs*100, "% and corresponds to", length(shuffle), "peptides.\n")
  
  if(length(shuffle) >= 2) {
  
    permutate <- c(shuffle[2:length(shuffle)], shuffle[1])
  
    MSRun[shuffle, parameters$QuantColnames] <- MSRun[permutate, parameters$QuantColnames]
  
    #Add annotation column for true and false positive identifications.
    MSRun$WrongID <- F
    MSRun$WrongID[shuffle] <- T
    
  } else {
    
    #Add annotation column for jsut true positive identifications.
    MSRun$WrongID <- F
    
  }
  
  cat("  - FDR addition finished.\n\n")
  
  # False PTM localization for different PTM types.
  MSRun$IsMisLocated <- F

  if (parameters$FracModProt > 0) {

    cat(" + PTM mis-localization:\n")
    
    for (mod in 1:length(parameters$PTMTypes)) {

      isModified <- which(sapply(MSRun$PTMType, function(x) any(x == parameters$PTMTypes[mod], na.rm=T)))
      modified <- MSRun[isModified, c("Sequence", "PTMPos", "PTMType")]

      # Count the number of modifiable residues
      tmpModPosCount <- stringr::str_locate_all(modified$Sequence, paste0(parameters$ModifiableResidues[[mod]], collapse="|"))
      ModifiableCount <- sapply(tmpModPosCount, function(x) length(x[,1]))
      ModCount <- sapply(modified$PTMType, function(x) sum(x == parameters$PTMTypes[mod]))

      # Positions in original table where modified peptide can get mislocated PTM.
      CanMisLoc <- sum(ModifiableCount > ModCount)

      # Number of modified peptides where we change localization.
      NumForMisLoc <- round(parameters$WrongLocalizations*CanMisLoc)

      # Mislocate PTMs for the modified peptides of NumForMisloc.
      if (NumForMisLoc > 0) {
        
        cat("  - For", parameters$PTMTypes[mod], "modification type,", NumForMisLoc, "modified peptides is selected for PTM re-location.\n")

        for (ind in sample(isModified[ModifiableCount > ModCount], NumForMisLoc)) {

          curr_pep <- MSRun[ind,]
          available_pos <- stringr::str_locate_all(curr_pep$Sequence, paste0(parameters$ModifiableResidues[[mod]],collapse="|"))[[1]][,1]
          PTMpos <- curr_pep$PTMPos[[1]][curr_pep$PTMType[[1]] == parameters$PTMTypes[mod]]

          remaining.pos <- available_pos[!(available_pos %in% PTMpos)]

          if(length(remaining.pos) != 0){

            # Dynamic number of mislocated PTMs: sample randomly from the actual PTMPos for a random size sampled from the length of the vector PTMpos.
            if(length(PTMpos) != 1){

              SubPTMpos <- sample(PTMpos, size = sample(1:length(PTMpos), size = 1, replace = F), replace = F)

            } else { # Unless, there is only a single PTM. Why not to use sample from a single integer: because sample function when data is a single number n creates a vector 1:n.

              SubPTMpos <- PTMpos

            }

            # From the remaining.pos sample with replacement a random number of PTMS with size of length SubPTMpos and keep the unique new positions.
            # Unique and sample by replace is used to cover the case when available positions are less than the actual PTM positions.
            if(length(remaining.pos) != 1){

              newPos <- unique(sample(remaining.pos, size = length(SubPTMpos), replace = T))

            } else {# Unless, there is only a single possible position.

              newPos <- remaining.pos

            }

            tPTMvec <- curr_pep$PTMPos[[1]]
            tPTMvec[sapply(SubPTMpos[1:length(newPos)], function(x) which(tPTMvec == x))] <- newPos
            curr_pep$PTMPos <- list(tPTMvec)
            curr_pep$IsMisLocated <- T
            MSRun[ind,] <- curr_pep

          }

        }
      }
    }
  }
  
  cat("  - PTM re-location finished.\n\n")
  
  if(!is.na(parameters$MaxNAPerPep) & parameters$MaxNAPerPep <= length(parameters$QuantColnames)){
    
    cat(" + Missing value filtering:\n")
    
    NAs <- apply(MSRun[, parameters$QuantColnames], 1, function(x) sum(is.na(x)))
    Which <- which(NAs <= parameters$MaxNAPerPep)
    
    cat("  - A total of", length(which), "peptides have removed, which have missing intensities in more than", parameters$MaxNAPerPep ,"samples.\n\n")
  
    if(length(which) > 0){
      
      MSRun <- MSRun[NAs <= parameters$MaxNAPerPep,]
      
    }
    
  }
  
  cat("#MS RUN SIMULATION - Finish\n\n")

  return(MSRun)
  
}
#####################

