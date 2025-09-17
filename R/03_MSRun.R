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
MSRunSim <- function(Digested, parameters, searchIndex = NULL) {

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
    # Shuffle the intensities of randomly selected peptides, to express the wrong identification.
    shuffle <- order(runif(nrow(MSRun)), decreasing = T)[seq_len(parameters$WrongIDs*nrow(MSRun))]
    message(" + Addition of false identification:")
    message("  - FDR selected is ", parameters$WrongIDs*100, "% and corresponds to ", length(shuffle), " peptides.")
    if (!is.null(searchIndex) && length(shuffle) > 0) {
        # Estimate donor windows and peptidoform space only if needed
        non_null <- vapply(searchIndex$proteins, function(x) !is.null(x), logical(1))
        nz_proteins <- sum(non_null)
        total_windows <- if (nz_proteins > 0) sum(vapply(searchIndex$proteins[non_null], function(x) x$n_valid, integer(1))) else 0L
        message("  - Donor windows (non-zero): ", total_windows, " across ", nz_proteins, " proteins.")
        # Estimate the size of the global peptidoform space (log-scale to avoid overflow)
        tot_log_pf <- .pm_total_peptidoforms_log(searchIndex, parameters)
        if (is.finite(tot_log_pf)) {
            tot_log10 <- tot_log_pf / log(10)
            expo <- floor(tot_log10)
            mant <- 10^(tot_log10 - expo)
            message(sprintf("  - Donor pool size (approx): %.3g x 10^%d peptidoforms.", mant, expo))
        }
    }

    # Build peptidoform donors from the global index if provided
    if (!is.null(searchIndex) && length(shuffle) > 0) {
        donors <- .pm_sample_uniform_peptidoforms(searchIndex, parameters, length(shuffle))
        if (!is.null(donors) && nrow(donors) == length(shuffle)) {
            # Preview a few substituted peptidoforms
            prev_n <- min(5L, nrow(donors))
            if (prev_n > 0) {
              ann <- mapply(.pm_annotate_peptidoform, donors$Peptide[seq_len(prev_n)], donors$PTMPos[seq_len(prev_n)], donors$PTMType[seq_len(prev_n)], SIMPLIFY = TRUE)
              message("  - Examples of substituted peptidoforms (Accession -> Peptidoform):")
              for (ii in seq_len(prev_n)) {
                message("    ", donors$Accession[[ii]], " -> ", ann[[ii]])
              }
            }
            # Keep intensities from MSRun; replace identities to simulate wrong IDs
            # Update sequence (I->L normalized) and PTM annotations; accession as list
            MSRun$Sequence[shuffle] <- gsub("[I]", "L", donors$Peptide)
            MSRun$PTMPos[shuffle] <- donors$PTMPos
            MSRun$PTMType[shuffle] <- donors$PTMType
            if ("Accession" %in% names(MSRun)) MSRun$Accession[shuffle] <- as.list(donors$Accession)
        }
    }

    # Annotate wrong IDs; keep existing intensities
    MSRun$WrongID <- FALSE
    if (length(shuffle) > 0) MSRun$WrongID[shuffle] <- TRUE

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

# Internal helpers to sample uniform peptidoforms from a full index
.pm_get_ptm_config <- function(parameters) {
  # Defaults if parameters are missing/NA
  default_types <- c("ph", "ox", "ac", "me", "de")
  default_residues <- list(
    ph = c("S", "T", "Y"),
    ox = c("M", "W", "C", "Y"),
    ac = c("K"),
    me = c("K", "R"),
    de = c("D", "E")
  )
  ptm_types <- tryCatch(parameters$PTMTypes, error = function(e) NULL)
  if (is.null(ptm_types) || length(ptm_types) == 0 || all(is.na(ptm_types))) {
    ptm_types <- list(default_types)
  }
  # Flatten one level if provided as list-of-one
  if (is.list(ptm_types) && length(ptm_types) == 1) ptm_types <- ptm_types[[1]]
  ptm_types <- ptm_types[!is.na(ptm_types)]
  modres <- tryCatch(parameters$ModifiableResidues, error = function(e) NULL)
  if (is.null(modres) || length(modres) == 0 || all(is.na(modres))) {
    modres_map <- default_residues
  } else {
    # If nested list, pick first layer
    if (is.list(modres) && length(modres) == 1 && is.list(modres[[1]])) {
      modres_map <- modres[[1]]
    } else if (is.list(modres)) {
      modres_map <- modres
    } else {
      modres_map <- default_residues
    }
  }
  list(types = ptm_types, residues = modres_map)
}

.pm_build_aa_maps <- function(parameters) {
  cfg <- .pm_get_ptm_config(parameters)
  aa_to_types <- setNames(vector("list", length = 26L), LETTERS)
  for (ptm in cfg$types) {
    aa <- cfg$residues[[ptm]]
    if (!is.null(aa) && length(aa) > 0) for (a in aa) aa_to_types[[a]] <- unique(c(aa_to_types[[a]], ptm))
  }
  aa_to_count <- setNames(integer(26L), LETTERS)
  for (a in names(aa_to_types)) aa_to_count[[a]] <- length(aa_to_types[[a]])
  list(aa_to_types = aa_to_types, aa_to_count = aa_to_count)
}

.pm_window_log_count <- function(seqi, start_pos, stop_pos, aa_to_count) {
  pep <- substring(seqi, start_pos, stop_pos)
  aa <- strsplit(pep, "", fixed = TRUE)[[1]]
  sum(log1p(as.numeric(aa_to_count[aa])))
}

.pm_sample_uniform_peptidoform_in_window <- function(seqi, start_pos, stop_pos, aa_to_types) {
  pep <- substring(seqi, start_pos, stop_pos)
  aa <- strsplit(pep, "", fixed = TRUE)[[1]]
  ptm_positions <- integer(0)
  ptm_types <- character(0)
  for (j in seq_along(aa)) {
    allowed <- aa_to_types[[aa[j]]]
    k <- length(allowed)
    draw <- sample.int(k + 1L, 1L) - 1L
    if (draw > 0L) {
      ptm_positions <- c(ptm_positions, j)
      ptm_types <- c(ptm_types, allowed[[draw]])
    }
  }
  list(
    Peptide = pep,
    Start = start_pos,
    Stop = stop_pos,
    PTMPos = if (length(ptm_positions)) list(ptm_positions) else list(integer(0)),
    PTMType = if (length(ptm_types)) list(ptm_types) else list(character(0))
  )
}

.pm_sample_uniform_peptidoforms <- function(index, parameters, N) {
  if (is.null(index) || is.na(N) || N <= 0) return(NULL)
  if (is.null(index$windows_flat) || is.null(index$window_probs) || length(index$window_probs) == 0) return(NULL)
  pick <- sample.int(nrow(index$windows_flat), size = N, replace = TRUE, prob = index$window_probs)
  out <- vector("list", length(pick))
  for (i in seq_along(pick)) {
    w <- index$windows_flat[pick[[i]], ]
    e <- index$proteins[[ w$p ]]
    sp <- w$start; ep <- w$stop
    spf <- .pm_sample_uniform_peptidoform_in_window(e$sequence, sp, ep, index$aa_to_types)
    out[[i]] <- data.frame(
      Accession = e$accession,
      Peptide = spf$Peptide,
      Start = spf$Start,
      Stop = spf$Stop,
      MC = w$mc,
      stringsAsFactors = FALSE
    )
    out[[i]]$PTMPos <- spf$PTMPos
    out[[i]]$PTMType <- spf$PTMType
  }
  do.call(rbind, out)
}

# Render a human-readable peptidoform by annotating PTMs in the sequence, e.g., S[ph]
.pm_annotate_peptidoform <- function(pep, pos, typ) {
  if (length(pos) == 0 || length(unlist(pos)) == 0) return(pep)
  p <- unlist(pos)
  t <- unlist(typ)
  chars <- strsplit(pep, "", fixed = TRUE)[[1]]
  ord <- order(p, decreasing = TRUE)
  for (i in ord) {
    k <- p[[i]]
    if (k >= 1 && k <= length(chars)) chars[k] <- paste0(chars[k], "[", t[[i]], "]")
  }
  paste0(chars, collapse = "")
}

# Compute total number of peptidoforms in log-scale (natural log)
.pm_total_peptidoforms_log <- function(index, parameters) {
  if (!is.null(index$total_log_pf)) return(index$total_log_pf)
  -Inf
}
#####################
