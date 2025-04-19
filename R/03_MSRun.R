################################################################################
#                               IN SILICO MS RUN                               #
################################################################################

#####################
#' Simulate MS run analysis
#'
#' This function simulates MS run analysis, which includes:
#' - Adding random noise to each sample due to MS instrument.
#' - Removing a specific proportion of peptides due to detection limitations.
#' - Removing a percentage of non-missing intensities based on probability weights depending on intensity values.
#' - Introducing peptide identification false discovery rate.
#' - Adding PTMs false localization rate for multiple PTM types.
#' - Filtering out peptides based on a maximum number of missing values threshold.
#' @param Digested A data frame containing digested peptides with associated data.
#' @param parameters A list containing various parameters for the MS run simulation.
#'
#' @return A data frame similar to `Digested` but with modifications introduced by the MS run simulation.
#'
#' @importFrom parallel detectCores makeCluster stopCluster setDefaultCluster clusterEvalQ parLapply
#' @importFrom stringr str_locate_all
#' @importFrom crayon red
#'
#' @keywords internal
MSRunSim <- function(Digested, parameters) {

    # TODO: what about multiples from different fractions?

    message("\n#MS RUN SIMULATION - Start\n")
    message(" + Noise addition:")
    message("  - The MS noise standard deviation is ", parameters$MSNoise, ".")

    # Introducing random noise to MS analysis due to MS instrument.
    matnoise <- matrix(rnorm(n = nrow(Digested) * length(parameters$QuantColnames), mean = 0, sd = parameters$MSNoise),
                       ncol=length(parameters$QuantColnames))

    Digested[ ,parameters$QuantColnames] <- Digested[ ,parameters$QuantColnames] + matnoise

    message("  - Noise added to all samples!\n")
    message(" + Detection limits:")
    message("  - Keeping ", parameters$PercDetectability*100, "% peptides with high detectability score (PeptideRanger).")

    # Sample a percentage of random peptides to be removed.
    if(parameters$PercDetectability < 1) {
        # get score threshold for lower percentage
        RFThreshold <- quantile(Digested$Detectability, 1-parameters$PercDetectability)
        remove <- Digested$Detectability <= RFThreshold
        MSRun <- Digested[!remove, ]
        message("  - A total of ", sum(remove), " peptides is removed with predicted detectability lower than ", RFThreshold, ".\n")

    } else {

        MSRun <- Digested
        message("  - No peptides were removed.\n")

    }

    #  message("  - The percentage of remaining peptides is", parameters$PercDetectedPep*100, "%.")

    # # Sample a percentage of random peptides to be removed.
    # if(parameters$PercDetectedPep != 1) {
    #
    #   remove <- sample(1:nrow(Digested), size = (1-parameters$PercDetectedPep)*nrow(Digested))
    #   MSRun <- Digested[-remove, ]
    #   message("  - A total of", (1-parameters$PercDetectedPep)*nrow(Digested), "peptides is removed.\n")
    #
    # } else {
    #
    #   MSRun <- Digested
    #   message("  - No peptides were removed.\n")
    #
    # }

    # Sample a random percentage of intensities to be removed.
    allVals <- as.vector(unlist(MSRun[ ,parameters$QuantColnames]))
    message(" + Removing peptide intensities:")
    message("  - The percentage of remaining intensities is ", parameters$PercDetectedVal*100, "% and the probabilities of selection depends on intensity values by ", parameters$WeightDetectVal, ".")
    myprob <-  (rank(-allVals)/length(allVals)) ^ parameters$WeightDetectVal

    #The probability for existing NA values from the above will be 1.
    #Thus, replacing probabilities corresponding to existing NAs to 0.
    myprob[is.na(allVals)] <- 0

    if(parameters$WeightDetectVal == 0){

        message(crayon::red("WARNING: missing values at peptide level (due to in silico MS detection) are determined at random."))
        message(crayon::red("         Change the parameter \'WeightDetectVal\' to a value > 0 to add weight to lowest intensities."))

    }

    #This is the fastest way to sample with replacement.
    #Reference: https://doi.org/10.1016/j.ipl.2005.11.003
    # method does not work as numbers become 0 -> move to log-scale
    #remove <- order(runif(length(allVals)) ^ (1/myprob), decreasing = T)[1:((1-parameters$PercDetectedVal)*length(allVals))]
    remove <- NULL
    if (parameters$PercDetectedVal < 1) {
        remove <- order(1/myprob * log(runif(length(allVals))), decreasing = T)[seq_len((1-parameters$PercDetectedVal)*length(allVals))]
        allVals[remove] <- NA
    }
    if (length(remove) > 0) {
        message("  - A total of ", length(remove), " intensities is removed.\n")
    } else {
        message("  - No intensities were removed.\n")
    }
    MSRun[ ,parameters$QuantColnames] <- matrix(allVals, ncol=length(parameters$QuantColnames))
    message("  - A total of ", length(remove), " intensities is removed.\n")
    # Shuffle the intensities of randomly selected peptides, to express the wrong identification.
    shuffle <- order(runif(nrow(MSRun)), decreasing = T)[seq_len(parameters$WrongIDs*nrow(MSRun))]
    message(" + Addition of false identification:")
    message("  - FDR selected is ", parameters$WrongIDs*100, "% and corresponds to ", length(shuffle), " peptides.")

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

    message("  - FDR addition finished.\n")

    # False PTM localization for different PTM types.
    MSRun$IsMisLocated <- F

    if (parameters$FracModProt > 0) {

        message(" + PTM mis-localization:")

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

                message("  - For ", parameters$PTMTypes[mod], " modification type, ", NumForMisLoc, " modified peptides are selected for PTM re-location.")

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

    message("  - PTM re-location finished.\n")

    if(!is.na(parameters$MaxNAPerPep) & parameters$MaxNAPerPep <= length(parameters$QuantColnames)){

        message(" + Missing value filtering:")

        NAs <- apply(MSRun[, parameters$QuantColnames], 1, function(x) sum(is.na(x)))
        Which <- which(NAs <= parameters$MaxNAPerPep)

        message("  - A total of ", length(which), " peptides have removed, which have missing intensities in more than ", parameters$MaxNAPerPep ," samples.\n")

        if(length(which) > 0){

            MSRun <- MSRun[NAs <= parameters$MaxNAPerPep,]

        }

    }

    message("#MS RUN SIMULATION - Finish\n")

    return(MSRun)

}
#####################

