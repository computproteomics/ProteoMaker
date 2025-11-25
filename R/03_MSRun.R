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
  # Require a valid search index; wrong-ID substitution depends on it
  if (is.null(searchIndex)) {
    stop("MSRunSim requires a non-NULL searchIndex. Build it with buildSearchIndexFromFasta() and pass it in.")
  }
  message(" + Noise addition:")
  message("  - The MS noise standard deviation is ", parameters$MSNoise, ".")
  # Report PTM sampling configuration for wrong-ID donors
  cfg <- .pm_get_ptm_config(parameters)
  if (length(cfg$types) == 0) {
    message("  - PTM sampling: none (donors are unmodified peptidoforms)")
  } else {
    message("  - PTM sampling types: ", paste(cfg$types, collapse = ", "))
  }

  # Introducing random noise to MS analysis due to MS instrument.
  matnoise <- matrix(rnorm(n = nrow(Digested) * length(parameters$QuantColnames), mean = 0, sd = parameters$MSNoise),
    ncol = length(parameters$QuantColnames)
  )

  Digested[, parameters$QuantColnames] <- Digested[, parameters$QuantColnames] + matnoise

  message("  - Noise added to all samples!\n")
  message(" + Detection limits:")
  message("  - Keeping ", parameters$PercDetectability * 100, "% peptides with high detectability score
          (PeptideRanger).")

  # Sample a percentage of random peptides to be removed.
  if (parameters$PercDetectability < 1) {
    # get score threshold for lower percentage
    RFThreshold <- quantile(Digested$Detectability, 1 - parameters$PercDetectability)
    remove <- Digested$Detectability <= RFThreshold
    MSRun <- Digested[!remove, ]
    message("  - A total of ", sum(remove), " peptides is removed with predicted detectability lower than ",
            RFThreshold, ".\n")
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
  allVals <- as.vector(unlist(MSRun[, parameters$QuantColnames]))
  message(" + Removing peptide intensities:")
  message("  - The percentage of remaining intensities is ", parameters$PercDetectedVal * 100,
          "% and the probabilities of selection depends on intensity values by ", parameters$WeightDetectVal, ".")
  myprob <- (rank(-allVals) / length(allVals))^parameters$WeightDetectVal

  # The probability for existing NA values from the above will be 1.
  # Thus, replacing probabilities corresponding to existing NAs to 0.
  myprob[is.na(allVals)] <- 0

  if (parameters$WeightDetectVal == 0) {
    message(crayon::red("WARNING: missing values at peptide level (due to in silico MS detection)
                        are determined at random."))
    message(crayon::red("         Change the parameter \'WeightDetectVal\' to a value > 0 to add weight to
                        lowest intensities."))
  }

  # This is the fastest way to sample with replacement.
  # Reference: https://doi.org/10.1016/j.ipl.2005.11.003
  # method does not work as numbers become 0 -> move to log-scale
  # remove <- order(runif(length(allVals)) ^ (1/myprob), decreasing = T)[1:((1-parameters$PercDetectedVal)*length(allVals))]
  remove <- NULL
  if (parameters$PercDetectedVal < 1) {
    remove <- order(1 / myprob * log(runif(length(allVals))), decreasing = T)[seq_len((1 - parameters$PercDetectedVal) *
                                                                                        length(allVals))]
    allVals[remove] <- NA
  }
  if (length(remove) > 0) {
    message("  - A total of ", length(remove), " intensities is removed.\n")
  } else {
    message("  - No intensities were removed.\n")
  }
  MSRun[, parameters$QuantColnames] <- matrix(allVals, ncol = length(parameters$QuantColnames))
  # Shuffle the intensities of randomly selected peptides, to express the wrong identification.
  shuffle <- order(runif(nrow(MSRun)), decreasing = T)[seq_len(parameters$WrongIDs * nrow(MSRun))]
  message(" + Addition of false identification:")
  message("  - FDR selected is ", parameters$WrongIDs * 100, "% and corresponds to ", length(shuffle), " peptides.")
  if (length(shuffle) > 0) {
    # Estimate donor windows and peptidoform space only if needed
    non_null <- vapply(searchIndex$proteins, function(x) !is.null(x), logical(1))
    nz_proteins <- sum(non_null)
    total_windows <- if (nz_proteins > 0) sum(vapply(searchIndex$proteins[non_null],
                                                     function(x) x$n_valid, integer(1))) else 0L
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
  if (length(shuffle) > 0) {
    donors <- .pm_sample_uniform_peptidoforms(searchIndex, parameters, length(shuffle))
    if (!is.null(donors) && nrow(donors) == length(shuffle)) {
      # Preview a few substituted peptidoforms
      prev_n <- min(5L, nrow(donors))
      if (prev_n > 0) {
        ann <- mapply(.pm_annotate_peptidoform, donors$Peptide[seq_len(prev_n)], donors$PTMPos[seq_len(prev_n)],
                      donors$PTMType[seq_len(prev_n)], SIMPLIFY = TRUE)
        message("  - Examples of substituted peptidoforms (Accession -> Peptidoform):")
        for (ii in seq_len(prev_n)) {
          message("    ", donors$Accession[[ii]], " -> ", ann[[ii]])
        }
      }
      # Diagnostics: how diverse are donor windows and sequences
      win_key <- paste(donors$Accession, donors$Start, donors$Stop, donors$MC, sep = ":")
      n_unique_windows <- length(unique(win_key))
      n_dupe_windows <- sum(duplicated(win_key))
      seq_norm <- gsub("[I]", "L", donors$Peptide)
      n_seq_dupes <- sum(duplicated(seq_norm))
      message("  - Donor diversity: ", n_unique_windows, " unique windows out of ", length(win_key),
              "; duplicate windows: ", n_dupe_windows, "; duplicate stripped sequences across donors: ",
              n_seq_dupes, ".")
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

  # Remove duplicate peptide sequences introduced by wrong-ID substitution
  if ("Sequence" %in% names(MSRun)) {
    dups <- duplicated(MSRun$Sequence)
    if (any(dups)) {
      removed <- sum(dups)
      MSRun <- MSRun[!dups, , drop = FALSE]
      message("  - Removed ", removed, " duplicate peptide sequences after wrong-ID step.")
    }
  }

  message("  - FDR addition finished.\n")

  # False PTM localization for different PTM types.
  MSRun$IsMisLocated <- F

  if (parameters$FracModProt > 0) {
    message(" + PTM mis-localization:")

    for (mod in 1:length(parameters$PTMTypes)) {
      isModified <- which(sapply(MSRun$PTMType, function(x) any(x == parameters$PTMTypes[mod], na.rm = T)))
      modified <- MSRun[isModified, c("Sequence", "PTMPos", "PTMType")]

      # Count the number of modifiable residues
      tmpModPosCount <- stringr::str_locate_all(modified$Sequence, paste0(parameters$ModifiableResidues[[mod]],
                                                                          collapse = "|"))
      ModifiableCount <- sapply(tmpModPosCount, function(x) length(x[, 1]))
      ModCount <- sapply(modified$PTMType, function(x) sum(x == parameters$PTMTypes[mod]))

      # Positions in original table where modified peptide can get mislocated PTM.
      CanMisLoc <- sum(ModifiableCount > ModCount)

      # Number of modified peptides where we change localization.
      NumForMisLoc <- round(parameters$WrongLocalizations * CanMisLoc)

      # Mislocate PTMs for the modified peptides of NumForMisloc.
      if (NumForMisLoc > 0) {
        message("  - For ", parameters$PTMTypes[mod], " modification type, ", NumForMisLoc,
                " modified peptides are selected for PTM re-location.")

        for (ind in sample(isModified[ModifiableCount > ModCount], NumForMisLoc)) {
          curr_pep <- MSRun[ind, ]
          available_pos <- stringr::str_locate_all(curr_pep$Sequence, paste0(parameters$ModifiableResidues[[mod]],
                                                                             collapse = "|"))[[1]][, 1]
          PTMpos <- curr_pep$PTMPos[[1]][curr_pep$PTMType[[1]] == parameters$PTMTypes[mod]]

          remaining.pos <- available_pos[!(available_pos %in% PTMpos)]

          if (length(remaining.pos) != 0) {
            # Dynamic number of mislocated PTMs: sample randomly from the actual PTMPos for a random size sampled from
            # the length of the vector PTMpos.
            if (length(PTMpos) != 1) {
              SubPTMpos <- sample(PTMpos, size = sample(1:length(PTMpos), size = 1, replace = F), replace = F)
            } else { # Unless, there is only a single PTM. Why not to use sample from a single integer: because sample
              # function when data is a single number n creates a vector 1:n.

              SubPTMpos <- PTMpos
            }

            # From the remaining.pos sample with replacement a random number of PTMS with size of length SubPTMpos and
            # keep the unique new positions.
            # Unique and sample by replace is used to cover the case when available positions are less than the actual
            # PTM positions.
            if (length(remaining.pos) != 1) {
              newPos <- unique(sample(remaining.pos, size = length(SubPTMpos), replace = T))
            } else { # Unless, there is only a single possible position.

              newPos <- remaining.pos
            }

            tPTMvec <- curr_pep$PTMPos[[1]]
            tPTMvec[sapply(SubPTMpos[1:length(newPos)], function(x) which(tPTMvec == x))] <- newPos
            curr_pep$PTMPos <- list(tPTMvec)
            curr_pep$IsMisLocated <- T
            MSRun[ind, ] <- curr_pep
          }
        }
      }
    }
  }

  message("  - PTM re-location finished.\n")

  if (!is.na(parameters$MaxNAPerPep) & parameters$MaxNAPerPep <= length(parameters$QuantColnames)) {
    message(" + Missing value filtering:")

    NAs <- apply(MSRun[, parameters$QuantColnames], 1, function(x) sum(is.na(x)))
    Which <- which(NAs <= parameters$MaxNAPerPep)

    message("  - A total of ", length(which), " peptides have removed, which have missing intensities in more than ",
            parameters$MaxNAPerPep, " samples.\n")

    if (length(which) > 0) {
      MSRun <- MSRun[NAs <= parameters$MaxNAPerPep, ]
    }
  }

  message("#MS RUN SIMULATION - Finish\n")

  return(MSRun)
}

# Internal helpers to sample uniform peptidoforms from a full index
.pm_get_ptm_config <- function(parameters) {
  ptm_types <- tryCatch(parameters$PTMTypes, error = function(e) NULL)
  # If not provided or all NA, treat as no PTMs
  if (is.null(ptm_types) || length(ptm_types) == 0 || all(is.na(ptm_types))) ptm_types <- character(0)
  # Flatten one level if provided as list-of-one
  if (is.list(ptm_types) && length(ptm_types) == 1) ptm_types <- ptm_types[[1]]
  ptm_types <- ptm_types[!is.na(ptm_types)]
  modres <- tryCatch(parameters$ModifiableResidues, error = function(e) NULL)
  if (is.null(modres) || length(modres) == 0 || all(is.na(modres))) {
    modres_map <- list()
  } else {
    # If nested list, pick first layer
    if (is.list(modres) && length(modres) == 1 && is.list(modres[[1]])) {
      modres_map <- modres[[1]]
    } else if (is.list(modres)) {
      modres_map <- modres
    } else {
      modres_map <- list()
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
  if (is.null(index) || is.na(N) || N <= 0) {
    return(NULL)
  }

  # Prefer per-protein sampling if win_count/pf_total are available
  has_pp <- length(index$proteins) > 0 && !is.null(index$proteins[[1]]$win_count)

  # Build PTM assignment map once
  aa_to_types <- if (!is.null(index$aa_to_types)) index$aa_to_types else .pm_build_aa_maps(parameters)$aa_to_types

  if (isTRUE(has_pp)) {
    # Protein weights
    pf <- vapply(index$proteins, function(e) if (is.null(e)) 0 else as.numeric(e$pf_total), numeric(1))
    if (!any(pf > 0)) {
      return(NULL)
    }
    prot_idx <- sample.int(length(index$proteins), size = N, replace = TRUE, prob = pf)
    by_prot <- split(seq_len(N), prot_idx)
    out <- vector("list", N)
    ptr <- 1L
    for (pi_chr in names(by_prot)) {
      pi <- as.integer(pi_chr)
      e <- index$proteins[[pi]]
      req <- length(by_prot[[pi_chr]])
      # Use log-counts when available to avoid overflow; otherwise fall back safely
      log_wc <- tryCatch(e$win_log_count, error = function(x) NULL)
      if (is.null(log_wc)) {
        wc <- e$win_count
        if (is.null(wc) || length(wc) == 0 || !any(is.finite(wc))) next
        # Stabilize by working in log-space
        log_wc <- log(pmax(wc, .Machine$double.xmin))
      }
      # Softmax in log-space for stable probabilities
      m <- max(log_wc)
      probs <- exp(log_wc - m)
      s <- sum(probs)
      if (!is.finite(s) || s <= 0) next
      probs <- probs / s
      sel <- sample.int(length(log_wc), size = req, replace = TRUE, prob = probs)
      for (jj in seq_len(req)) {
        k <- sel[[jj]]
        sp <- e$win_start[[k]]
        ep <- e$win_stop[[k]]
        spf <- .pm_sample_uniform_peptidoform_in_window(e$sequence, sp, ep, aa_to_types)
        row <- data.frame(
          Accession = e$accession,
          Peptide = spf$Peptide,
          Start = spf$Start,
          Stop = spf$Stop,
          MC = e$win_mc[[k]],
          stringsAsFactors = FALSE
        )
        row$PTMPos <- spf$PTMPos
        row$PTMType <- spf$PTMType
        out[[ptr]] <- row
        ptr <- ptr + 1L
      }
    }
    return(do.call(rbind, out))
  }

  # No per-protein data available
  return(NULL)
}

# Render a human-readable peptidoform by annotating PTMs in the sequence, e.g., S[ph]
.pm_annotate_peptidoform <- function(pep, pos, typ) {
  if (length(pos) == 0 || length(unlist(pos)) == 0) {
    return(pep)
  }
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
  if (is.null(index) || is.null(index$proteins)) {
    return(-Inf)
  }
  tot <- sum(vapply(index$proteins,
                    function(e) if (is.null(e) || is.null(e$pf_total)) 0 else as.numeric(e$pf_total), numeric(1)))
  if (!is.finite(tot) || tot <= 0) {
    return(-Inf)
  }
  log(tot)
}
#####################
