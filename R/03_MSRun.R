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
MSRunSim <- function(Digested, parameters, search_index) {

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
    # Select rows to be assigned wrong identifications
    n_wrong <- as.integer(parameters$WrongIDs * nrow(MSRun))
    wrong_rows <- integer(0)
    if (n_wrong > 0) {
      wrong_rows <- order(runif(nrow(MSRun)), decreasing = TRUE)[seq_len(n_wrong)]
    }
    message(" + Addition of false identification:")
    message("  - FDR selected is ", parameters$WrongIDs*100, "% and corresponds to ", length(wrong_rows), " peptides.")

    # Initialize WrongID flag
    MSRun$WrongID <- FALSE

    if (length(wrong_rows) > 0) {
      # Prepare index sampling and candidate pool
      nprot <- length(search_index$proteins)
      if (nprot <= 0) stop("No decoy proteins available in index for wrong-ID assignment.")
      w <- search_index$weights
      if (length(w) != nprot || any(!is.finite(w)) || sum(w, na.rm = TRUE) <= 0) {
        w <- rep(1 / nprot, nprot)
      }
      tmpl_cache <- vector("list", nprot)
      ok_cache <- vector("list", nprot)
      # Maintain normalized sequences to avoid duplicates
      existing_seq_norm <- as.character(gsub("[I]", "L", MSRun$Sequence))
      newly_assigned <- character(0)
      unique_pool <- new.env(parent = emptyenv())
      # Build per-protein candidate lists; dismiss problematic proteins
      for (pj in seq_len(nprot)) {
        ip <- search_index$proteins[[pj]]
        tmpl <- NULL
        if (!is.null(ip$starts) && !is.null(ip$stops) && !is.null(ip$valid_mc)) {
          tmpl <- protein_peptide_template_from_index(ip)
        }
        if (!is.null(tmpl) && nrow(tmpl) > 0) {
          cn <- gsub("[I]", "L", tmpl$Peptide)
          ok <- which(!(cn %in% existing_seq_norm))
          if (length(ok) > 0) {
            tmpl_cache[[pj]] <- tmpl
            ok_cache[[pj]] <- ok
            for (uu in unique(cn[ok])) assign(uu, TRUE, envir = unique_pool)
          } else {
            ok_cache[[pj]] <- integer(0)
          }
        } else {
          ok_cache[[pj]] <- integer(0)
        }
      }
      has_cand <- vapply(ok_cache, function(x) length(x) > 0, logical(1))
      if (!any(has_cand)) stop("No decoy peptides available for wrong-ID assignment (all candidates overlap current dataset).")
      w[!has_cand] <- 0
      if (sum(w) > 0) w <- w / sum(w) else stop("No decoy proteins with candidates available after filtering.")
      total_unique_decoys <- length(ls(envir = unique_pool, all.names = TRUE))
      message("  - Wrong-ID decoy pool: ", total_unique_decoys, " unique peptides across ", sum(has_cand), " proteins.")

      for (ri in wrong_rows) {
        # Try sampling until a protein with remaining candidates is found
        attempts <- 0
        repeat {
          attempts <- attempts + 1
          if (sum(w) == 0) stop("No decoy peptides remaining for wrong-ID assignment (assigned ", length(newly_assigned), " of ", length(wrong_rows), ").")
          pick_idx <- sample.int(nprot, size = 1, replace = TRUE, prob = w)
          if (is.na(pick_idx) || pick_idx < 1 || pick_idx > nprot) next
          tmpl <- tmpl_cache[[pick_idx]]
          if (is.null(tmpl) || nrow(tmpl) == 0) { w[pick_idx] <- 0; if (sum(w)>0) w <- w / sum(w); next }
          cn_all <- gsub("[I]", "L", tmpl$Peptide)
          ok <- ok_cache[[pick_idx]]
          if (length(ok) == 0) { w[pick_idx] <- 0; if (sum(w)>0) w <- w / sum(w); next }
          ok2 <- ok[!(cn_all[ok] %in% newly_assigned)]
          if (length(ok2) == 0) { w[pick_idx] <- 0; if (sum(w)>0) w <- w / sum(w); next }
          pick_idx2 <- sample(ok2, size = 1)
          decoy <- tmpl[pick_idx2, ]
          decoy_seq_norm <- cn_all[pick_idx2]
          # consume this candidate
          ok_cache[[pick_idx]] <- setdiff(ok, pick_idx2)
          if (length(ok_cache[[pick_idx]]) == 0) { w[pick_idx] <- 0; if (sum(w)>0) w <- w / sum(w) }
          break
        }
        # Overwrite identification fields; keep intensities unchanged
        MSRun$Sequence[ri] <- decoy_seq_norm
        MSRun$Peptide[ri] <- list(decoy$Peptide)
        MSRun$Start[ri] <- list(decoy$Start)
        MSRun$Stop[ri]  <- list(decoy$Stop)
        MSRun$MC[ri]    <- list(decoy$MC)
        MSRun$MZ1[ri]   <- decoy$MZ1
        MSRun$MZ2[ri]   <- decoy$MZ2
        MSRun$MZ3[ri]   <- decoy$MZ3
        # Assign precise peptide->protein mapping if available in index
        if (!is.null(search_index$pep2prot) && length(search_index$pep2prot) > 0 && decoy_seq_norm %in% names(search_index$pep2prot)) {
          MSRun$Accession[ri] <- list(search_index$pep2prot[[decoy_seq_norm]])
        } else {
          MSRun$Accession[ri] <- list(idx_prot$accession)
        }
        # Assign PTMs on decoy if the original identification carried PTMs and decoy has compatible residues
        orig_ptm_types <- MSRun$PTMType[[ri]]
        if (!is.null(orig_ptm_types) && length(orig_ptm_types) > 0) {
          # Determine available modifiable positions on the decoy segment for each PTM type
          decoy_start <- decoy$Start
          decoy_stop  <- decoy$Stop
          assigned_pos <- integer(0)
          assigned_type <- character(0)

          # Mass shifts per PTM type (if available)
          ptm_mass_map <- NULL
          if (!is.null(parameters$PTMTypesMass)) {
            ptm_mass_map <- parameters$PTMTypesMass[[1]]
            if (!is.null(parameters$PTMTypes) && length(parameters$PTMTypes) > 0) {
              names(ptm_mass_map) <- parameters$PTMTypes[[1]]
            }
          }

          for (t in unique(orig_ptm_types)) {
            # how many of this type were on the original
            n_needed <- sum(orig_ptm_types == t)
            avail_sites <- integer(0)
            # Prefer index-provided protein-level modifiable sites per PTM type
            if (!is.null(idx_prot$modpos) && !is.null(idx_prot$modpos[[t]])) {
              avail_sites <- idx_prot$modpos[[t]]
              # restrict to this peptide window
              if (length(avail_sites) > 0) {
                avail_sites <- avail_sites[avail_sites >= decoy_start & avail_sites <= decoy_stop]
              }
            }
            # Fallback: scan peptide sequence for modifiable residues by type
            if (length(avail_sites) == 0 && !is.null(parameters$ModifiableResidues)) {
              mod_aas <- NULL
              # handle both list-of-list structure and flat mapping gracefully
              if (!is.null(parameters$ModifiableResidues[[1]]) && !is.null(parameters$ModifiableResidues[[1]][[t]])) {
                mod_aas <- parameters$ModifiableResidues[[1]][[t]]
              } else if (!is.null(parameters$ModifiableResidues[[t]])) {
                mod_aas <- parameters$ModifiableResidues[[t]]
              }
              if (!is.null(mod_aas)) {
                seg <- substring(idx_prot$sequence, decoy_start, decoy_stop)
                if (nchar(seg) > 0) {
                  seg_chars <- strsplit(seg, "")[[1]]
                  rel <- which(seg_chars %in% mod_aas)
                  if (length(rel) > 0) {
                    avail_sites <- decoy_start + rel - 1
                  }
                }
              }
            }
            if (length(avail_sites) > 0) {
              n_pick <- min(n_needed, length(avail_sites))
              pick <- if (n_pick > 0) sample(avail_sites, size = n_pick, replace = FALSE) else integer(0)
              if (length(pick) > 0) {
                assigned_pos <- c(assigned_pos, pick - decoy_start + 1)
                assigned_type <- c(assigned_type, rep(t, length(pick)))
              }
            }
          }

          if (length(assigned_pos) > 0) {
            # Assign PTM positions/types to decoy
            MSRun$PTMPos[ri]  <- list(as.integer(assigned_pos))
            MSRun$PTMType[ri] <- list(assigned_type)
            # Adjust M/Z by the total mass shift from assigned PTMs, if mass map is available
            if (!is.null(ptm_mass_map)) {
              mass_add <- sum(vapply(assigned_type, function(tt) {
                mv <- ptm_mass_map[[tt]]
                if (is.null(mv) || is.na(mv)) 0 else as.numeric(mv)
              }, numeric(1)))
              if (!is.na(mass_add) && mass_add != 0) {
                MSRun$MZ1[ri] <- MSRun$MZ1[ri] + mass_add
                MSRun$MZ2[ri] <- MSRun$MZ2[ri] + mass_add/2
                MSRun$MZ3[ri] <- MSRun$MZ3[ri] + mass_add/3
              }
            }
          } else {
            # no compatible sites on decoy; leave as unmodified
            MSRun$PTMPos[ri]  <- list(NULL)
            MSRun$PTMType[ri] <- list(NULL)
          }
        } else {
          # Original had no PTMs; keep decoy unmodified
          MSRun$PTMPos[ri]  <- list(NULL)
          MSRun$PTMType[ri] <- list(NULL)
        }
        # Flag as wrong ID
        MSRun$WrongID[ri] <- TRUE

        # Track assigned decoy to avoid duplicates for subsequent rows
        newly_assigned <- c(newly_assigned, decoy_seq_norm)
        existing_seq_norm <- c(existing_seq_norm, decoy_seq_norm)
      }
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
