################################################################################
#               Data analysis of the outputs of the peptide quan.              #
################################################################################

#' Summarize Protein-Group Abundance from Peptidoform-Level Data
#'
#' This function summarizes protein-group abundance from peptidoform-level data.
#' It optionally removes modified peptidoforms, groups peptidoforms by their
#' parent protein groups, and applies a summarization method (e.g., sum of top 3,
#' median polish) to estimate protein-group abundance. The function supports
#' parallel processing to speed up the summarization process.
#'
#' @param peptable A data frame containing peptidoform-level data, including accession numbers,
#' sequence, PTM types, and quantification columns.
#' @param parameters A list of parameters, including:
#' \describe{
#'   \item{ProtSummarization}{The method used for protein summarization, e.g., "sum.top3" or "medpolish".}
#'   \item{MinUniquePep}{The minimum number of unique peptides required to summarize a protein.}
#'   \item{QuantColnames}{The names of the columns containing quantification data.}
#'   \item{Cores}{The number of cores to use for parallel processing.}
#'   \item{ClusterType}{The type of cluster to use for parallel processing (e.g., "FORK", "PSOCK").}
#' }
#'
#' @return A data frame containing summarized protein-group data, where each row represents a
#' protein group, and the columns include protein-group information and summarized quantification data.
#'
#' @importFrom parallel detectCores makeCluster setDefaultCluster clusterExport parLapply stopCluster
#' @importFrom tidyr pivot_longer
#' @importFrom MASS rlm
#'
#' @keywords internal
#'
proteinSummarisation <- function(peptable, parameters) {
  method <- parameters$ProtSummarization

  minUniquePep <- parameters$MinUniquePep

  includeModPep <- parameters$IncludeModPep

  sharedPep <- parameters$SharedPep

  QuantColnames <- parameters$QuantColnames

  if (!includeModPep) {
    # Remove all modified peptidoforms
    peptable <- peptable[!sapply(peptable$PTMType, function(x) length(x) > 0), ]
    message("  - Removed all modified peptidoforms, remaining number of peptidoforms: ", nrow(peptable), "")
  } else {
    message("  - Keeping modified peptidoforms")
  }
  message("  - Remaining number of peptidoforms: ", nrow(peptable), "")

  message(" + Protein summarisation")

  message("  - Minimum unique peptidoforms per protein group: ", minUniquePep, "")
  message("  - Protein-group summarisation using the ", method, " approach.")

  # writing new column with unlisted and merged protein names
  peptable$merged_accs <- sapply(peptable$Accession, function(x) paste(sort(unique(unlist(x))), collapse = ";"))
  peptable$num_accs <- sapply(peptable$Accession, function(x) length(unique(x)))

  # Sort table according to protein accession, needs to stay in this order!
  peptable <- peptable[order(peptable$merged_accs), ]
  message("  - Sorted protein table")

  # Reducing table to relevant columns
  if (!sharedPep) {
    peptable <- peptable[peptable$num_accs == 1, ]
  }

  pep_key <- if ("Peptidoform" %in% names(peptable)) peptable$Peptidoform else peptable$Sequence

  # Vector with row indices of protein groups
  all_accs <- peptable$merged_accs
  prot_ind <- 1
  names(prot_ind) <- all_accs[1]
  for (i in 2:nrow(peptable)) {
    if (all_accs[i - 1] != all_accs[i]) {
      prot_ind <- c(prot_ind, i)
      names(prot_ind)[length(prot_ind)] <- all_accs[i]
    }
  }
  prot_ind <- c(prot_ind, nrow(peptable))
  other_cols <- colnames(peptable)[!colnames(peptable) %in% QuantColnames]
  message("  - built protein index for faster summarization")

  # Diagnostics: duplicated stripped sequences within each protein group
  dup_group_info <- list()
  dup_total_rows <- 0L
  dup_group_count <- 0L
  n_groups <- length(prot_ind) - 1L
  if (n_groups > 0L) {
    for (gi in seq_len(n_groups)) {
      start_i <- prot_ind[gi]
      end_i <- prot_ind[gi + 1L] - 1L
      if (end_i >= start_i) {
        seqs <- pep_key[start_i:end_i]
        wrongid <- if ("WrongID" %in% names(peptable)) peptable$WrongID[start_i:end_i] else rep(NA, length(seqs))
        dups <- seqs[duplicated(seqs)]
        if (length(dups) > 0L) {
          dup_group_count <- dup_group_count + 1L
          dup_total_rows <- dup_total_rows + length(dups)
          u <- unique(dups)
          # total counts per duplicated sequence
          ct <- vapply(u, function(s) sum(seqs == s), numeric(1))
          # how many of those rows are WrongID
          wt <- vapply(u, function(s) sum(wrongid[seqs == s], na.rm = TRUE), numeric(1))
          # store as a small data.frame for later pretty printing
          dup_group_info[[names(prot_ind)[gi]]] <- data.frame(seq = u, n = as.integer(ct), wrongID = as.integer(wt), stringsAsFactors = FALSE)
        }
      }
    }
    if (dup_group_count > 0L) {
      message(
        "  - Duplicated stripped sequences within protein groups: ", dup_total_rows,
        " duplicate rows across ", dup_group_count, " groups."
      )
      max_groups <- 5L
      max_each <- 5L
      shown <- head(names(dup_group_info), max_groups)
      for (g in shown) {
        df <- dup_group_info[[g]]
        o <- order(df$n, decreasing = TRUE)
        df <- df[o, , drop = FALSE]
        if (nrow(df) > max_each) df <- df[seq_len(max_each), , drop = FALSE]
        parts <- paste0(df$seq, " (n=", df$n, ", wrongID=", df$wrongID, ")")
        message("    ", g, ": ", paste(parts, collapse = "; "))
      }
      if (length(dup_group_info) > max_groups) {
        message("    ... ", length(dup_group_info) - max_groups, " more groups show duplicates.")
      }
    } else {
      message("  - No duplicated stripped sequences within protein groups detected.")
    }
  }

  # Initiate and fill matrix with proteins
  protmat <- as.data.frame(matrix(ncol = ncol(peptable), nrow = length(prot_ind)))
  rownames(protmat) <- names(prot_ind)
  colnames(protmat) <- colnames(peptable)

  message("  - Initiated protein matrix for ", length(prot_ind), " protein groups")
  # Diagnostics: potential duplicate first sequences across protein groups
  if ((length(prot_ind) - 1) > 1 && length(pep_key) > 0) {
    first_rows <- prot_ind[1:(length(prot_ind) - 1)]
    first_seq <- pep_key[first_rows]
    dup_count <- sum(duplicated(first_seq))
    if (dup_count > 0) {
      message("  - Note: ", dup_count, " protein groups share the same first peptide sequence; row names will use group keys to avoid clashes.")
    }
  }
  message("  - Summarizing protein groups, this can take a while")

  # Function to summarize protein groups
  summarizeProtein <- function(tmp) {
    out <- NULL
    used_fallback <- FALSE
    if (nrow(tmp) >= minUniquePep) {
      tmp <- as.matrix(tmp)
      if (method == "sum.top3") {
        tmp <- tmp[order(rowSums(tmp), decreasing = T), , drop = F]
        if (nrow(tmp) >= 3) {
          out <- log2(colSums(2^tmp[1:3, ], na.rm = T))
        } else {
          out <- log2(colSums(2^tmp, na.rm = T))
        }
      } else if (method == "median") {
        out <- apply(tmp, 2, median, na.rm = T)
      } else if (method == "mean") {
        out <- apply(tmp, 2, mean, na.rm = T)
      } else if (method == "sum") {
        out <- log2(colSums(2^tmp, na.rm = T))
      } else if (method == "medpolish") {
        summed <- NULL
        if (nrow(tmp) <= 3) {
          summed <- apply(tmp, 2, median, na.rm = TRUE)
        } else {
          summed <- medpolish(tmp, na.rm = T, trace.iter = F)$col
        }
        if (length(summed) > 0) {
          out <- summed
        }
      } else if (method == "rlm") {
        if (nrow(tmp) > 1) {
          tmp <- as.data.frame(tmp)
          tmp$peptide <- rownames(tmp)
          long_df <- pivot_longer(tmp,
                                  cols = -peptide,
                                  names_to = "sample",
                                  values_to = "intensity"
          )
          out <- setNames(rep(NA_real_, length(QuantColnames)), QuantColnames)
          long_df <- long_df[!is.na(long_df$intensity), , drop = FALSE]
          if (nrow(long_df) > 0) {
            long_df$peptide <- factor(long_df$peptide)
            long_df$sample <- factor(long_df$sample, levels = QuantColnames)
            long_df$sample <- droplevels(long_df$sample)
            n_pep <- nlevels(long_df$peptide)
            n_samp <- nlevels(long_df$sample)
            n_param <- n_pep + n_samp - 1L
            if (n_pep > 1 && n_samp > 1 && nrow(long_df) > n_param) {
              mm <- model.matrix(~ peptide + sample, data = long_df)
              if (qr(mm)$rank == ncol(mm)) {
                fit <- rlm(intensity ~ peptide + sample, data = long_df, na.action = na.omit)
                coefs <- coef(fit)
                sample_coefs <- coefs[grep("^sample", names(coefs))]

                # Include intercept (baseline)
                baseline <- coefs["(Intercept)"]
                sample_names <- gsub("^sample", "", names(sample_coefs))
                baseline_sample <- levels(long_df$sample)[1]
                out[baseline_sample] <- baseline
                out[sample_names] <- baseline + sample_coefs
              }
            }
          }
          if (all(is.na(out))) {
            out <- apply(as.matrix(tmp[, QuantColnames, drop = FALSE]), 2, mean, na.rm = TRUE)
            used_fallback <- TRUE
          }
        } else {
          out <- tmp[1, QuantColnames]
        }
      } else {
        stop("No valid method provided!")
      }
    }
    attr(out, "used_fallback") <- used_fallback
    return(out)
  }

  if (!is.null(parameters$Cores)) {
    cores <- parameters$Cores

    if (parallel::detectCores() <= parameters$Cores) {
      cores <- parallel::detectCores() - 1
    }

    cluster <- parallel::makeCluster(cores, type = parameters$ClusterType)
    parallel::setDefaultCluster(cluster)
    parallel::clusterExport(cluster, c(
      "peptable", "summarizeProtein", "minUniquePep",
      "prot_ind", "other_cols", "QuantColnames"
    ),
    envir = environment()
    )
    # load required packages on all workers
    parallel::clusterEvalQ(cluster, {
      library(MASS)
      # Add others if needed
      library(dplyr)
      library(tidyr)
      # etc.
    })
    proteins <- parallel::parLapply(cluster, 1:(length(prot_ind) - 1), function(i) {
      tmp <- as.data.frame(peptable[prot_ind[i]:(prot_ind[i + 1] - 1), ])
      key <- if ("Peptidoform" %in% names(tmp)) tmp$Peptidoform else tmp$Sequence
      rownames(tmp) <- make.unique(key)
      out <- tmp[1, , drop=F]
      tout <- summarizeProtein(tmp[, QuantColnames, drop = F])
      used_fallback <- isTRUE(attr(tout, "used_fallback"))
      if (!is.null(tout)) {
        out[QuantColnames] <- tout
        # add other information
        out[other_cols] <- sapply(tmp[, other_cols], function(x) paste(unlist(x), collapse = ";"))
        # ensure unique row name per protein group (use group key)
        rownames(out) <- names(prot_ind)[i]
      } else {
        out <- NULL
      }
      return(list(out = out, used_fallback = used_fallback))
    })
    parallel::stopCluster(cluster)
  } else {
    proteins <- lapply(1:(length(prot_ind) - 1), function(i) {
      tmp <- as.data.frame(peptable[prot_ind[i]:(prot_ind[i + 1] - 1), ])
      key <- if ("Peptidoform" %in% names(tmp)) tmp$Peptidoform else tmp$Sequence
      rownames(tmp) <- make.unique(key)
      out <- tmp[1, ]
      tout <- summarizeProtein(tmp[, QuantColnames, drop = F])
      used_fallback <- isTRUE(attr(tout, "used_fallback"))
      if (!is.null(tout)) {
        out[QuantColnames] <- tout
        # add other information
        out[other_cols] <- sapply(tmp[, other_cols], function(x) paste(unlist(x), collapse = ";"))
        rownames(out) <- names(prot_ind)[i]
      } else {
        out <- NULL
      }
      return(list(out = out, used_fallback = used_fallback))
    })
  }
  fallback_count <- sum(vapply(proteins, function(x) isTRUE(x$used_fallback), logical(1)))
  proteins <- lapply(proteins, function(x) x$out)
  # join all protein data
  proteins <- Filter(Negate(is.null), proteins)
  protmat <- do.call(rbind, proteins)
  protmat[protmat == -Inf] <- NA
  protmat <- protmat[rowSums(is.na(protmat[, QuantColnames])) < length(QuantColnames), ]

  cat("  - Finished summarizing into ", nrow(protmat), " protein groups")
  if (method == "rlm" && fallback_count > 0) {
    cat("\n  - rlm fallback to mean used in ", fallback_count, " protein groups")
  }
  #  for (i in parameters$QuantColnames) protmat[,i] <- as.numeric(protmat[,i])

  return(protmat)
}


#' Calculate PTM Site Occupancy
#'
#' Estimates the per-condition PTM-site occupancy (stoichiometry) using a
#' mass-balance three-ratio approach.  For each modified peptidoform the
#' method uses three between-condition fold-change ratios:
#' \describe{
#'   \item{\eqn{R_{m,c}}}{Fold-change of the \emph{modified} peptidoform,
#'     \eqn{R_{m,c} = 2^{\bar{l}_{m,c} - \bar{l}_{m,1}}}.}
#'   \item{\eqn{R_{u,c}}}{Fold-change of the \emph{same-sequence unmodified}
#'     counterpart peptidoform,
#'     \eqn{R_{u,c} = 2^{\bar{l}_{u,c} - \bar{l}_{u,1}}}.}
#'   \item{\eqn{R_{prot,c}}}{Fold-change of the protein, estimated from the
#'     geometric mean (in linear space) of all unmodified peptides from the
#'     same protein accession whose sequence does \emph{not} appear in any
#'     modified form in the data,
#'     \eqn{R_{prot,c} = 2^{\bar{l}_{prot,c} - \bar{l}_{prot,1}}}.}
#' }
#' Here \eqn{\bar{l}_{x,c}} is the mean log2 intensity across replicates of
#' condition \eqn{c} and subscript 1 denotes the reference condition.  For the
#' protein ratio, log2 intensities from multiple non-counterpart unmodified rows
#' are pooled first via \code{colMeans} (geometric mean in linear space) and
#' then averaged across replicates.  Missing values are removed
#' (\code{na.rm = TRUE}) at every averaging step; consequently, an individual
#' missing replicate is silently excluded from the mean for that condition.
#' However, if \emph{all} replicates of a condition are missing, the condition
#' mean becomes \code{NaN}.
#'
#' Under a simple mass-balance model, the reference-condition occupancy satisfies
#' the identity \eqn{occ_1 = (R_{prot,c} - R_{u,c}) / (R_{m,c} - R_{u,c})} for
#' every non-reference condition \eqn{c}.  For robustness this quantity is
#' estimated as the mean over all \eqn{C - 1} non-reference conditions:
#'
#' \deqn{occ_1 = \frac{1}{C - 1} \sum_{c=2}^{C}
#'   \frac{R_{prot,c} - R_{u,c}}{R_{m,c} - R_{u,c}}}
#'
#' Per-condition occupancy is then obtained by:
#'
#' \deqn{occ_c = occ_1 \cdot \frac{R_{m,c}}{R_{prot,c}}, \quad c = 2, \ldots, C}
#'
#' The reference condition column (\code{C_1}) contains \eqn{occ_1}; all other
#' columns contain \eqn{occ_c}.  When \eqn{R_{m,c} = R_{u,c}} for any
#' condition (denominator of the \eqn{occ_1} estimator is zero), the
#' per-condition term yields \code{NaN}, which propagates through the mean to
#' \eqn{occ_1} and consequently to \emph{all} output columns.  A peptide
#' sequence is omitted from the result when: (a) no same-sequence unmodified
#' row is present; (b) no non-counterpart unmodified row from the same protein
#' is available to estimate \eqn{R_{prot}}; or (c) more than one modified or
#' more than one unmodified row exists for that sequence (a warning is issued
#' and the sequence is skipped). If peptide start positions are available,
#' modified peptides with more than one possible protein start position are
#' also skipped because their protein-level PTM position is ambiguous.
#'
#' @param peptable A data frame of peptidoform-level data as produced by the
#'   ProteoMaker pipeline (output of \code{MSRunSim} or \code{runPolySTest}).
#'   Required columns: \code{Sequence} (character), \code{PTMType} (list),
#'   \code{PTMPos} (list), \code{Accession} (list), and one column per sample
#'   as given by \code{parameters$QuantColnames}.  Quantification values must
#'   be on the log2 scale. If \code{Start} is present, it is used to report
#'   protein-level PTM positions in \code{ProteinPTMPos}.
#' @param parameters A named list of analysis parameters.  Must contain:
#'   \describe{
#'     \item{QuantColnames}{Character vector of column names holding per-sample
#'       log2 intensities.  Samples must be ordered as all replicates of
#'       condition 1 first, then all replicates of condition 2, etc.}
#'     \item{NumCond}{Integer.  Number of experimental conditions (\eqn{\geq 2}).}
#'     \item{NumReps}{Integer.  Number of replicates per condition.}
#'   }
#'
#' @return A data frame with one row per modified peptidoform that has both a
#'   same-sequence unmodified counterpart and at least one non-counterpart
#'   unmodified peptide from the same protein.  Columns are:
#'   \describe{
#'     \item{Sequence}{Stripped peptide sequence.}
#'     \item{Accession}{Protein accession(s) (list column).}
#'     \item{PTMPos}{Modification positions within the peptide (list column).}
#'     \item{ProteinPTMPos}{Modification positions within the protein sequence
#'       (list column), or \code{NA} when peptide start positions are not
#'       available.}
#'     \item{PTMType}{Modification types (list column).}
#'     \item{C_1, C_2, \ldots, C_NumCond}{One numeric column per condition.
#'       \code{C_1} contains the estimated reference-condition occupancy
#'       \eqn{occ_1}.  Columns \code{C_2} through \code{C_NumCond} contain
#'       \eqn{occ_c = occ_1 \cdot R_{m,c} / R_{prot,c}}.  Values may fall
#'       outside \eqn{[0, 1]} or be \code{NaN} when the data are inconsistent
#'       with the mass-balance model (e.g.\ \eqn{R_{m,c} = R_{u,c}}).}
#'   }
#'   Returns an empty \code{data.frame} when no modified peptides are present,
#'   when \code{NumCond < 2}, or when no qualifying modified peptide is found.
#'
#' @references
#' Sharma, K. et al. (2014) Ultradeep Human Phosphoproteome Reveals a Distinct
#' Regulatory Nature of Tyr and Ser/Thr-Based Signaling.
#' \emph{Cell Reports}, \bold{8}(5), 1583--1594.
#' \doi{10.1016/j.celrep.2014.07.036}
#'
#' Olsen, J.V. et al. (2010) Quantitative Phosphoproteomics Reveals Widespread
#' Full Phosphorylation Site Occupancy During Mitosis.
#' \emph{Science Signaling}, \bold{3}(104), ra3.
#' \doi{10.1126/scisignal.2000475}
#'
#' @examples
#' \dontrun{
#' # After running the ProteoMaker pipeline and obtaining StatsPep:
#' occ <- calcPTMOccupancy(StatsPep, Param)
#' head(occ)
#' }
#'
#' @importFrom mvtnorm pmvnorm
#'
#' @export
calcPTMOccupancy <- function(peptable, parameters) {
  QuantColnames <- parameters$QuantColnames
  NumCond       <- parameters$NumCond
  NumReps       <- parameters$NumReps

  has_ptm <- lengths(peptable$PTMType) > 0

  if (sum(has_ptm) == 0) {
    message("calcPTMOccupancy: no modified peptides found; returning empty table.")
    return(data.frame())
  }

  if (is.null(NumCond) || is.null(NumReps) || NumCond < 2L) {
    message("calcPTMOccupancy: NumCond >= 2 required for between-condition occupancy; returning empty table.")
    return(data.frame())
  }

  message(" + Calculating PTM site occupancy")

  unique_values <- function(x) unique(unlist(x))

  # Group QuantColnames into per-condition replicate blocks
  cond_cols      <- lapply(seq_len(NumCond), function(c) {
    QuantColnames[seq.int((c - 1L) * NumReps + 1L, c * NumReps)]
  })
  out_cond_names <- paste0("C_", seq_len(NumCond))
  out_comp_names <- paste0("prob_", out_cond_names[seq_len(NumCond-1)+1], "_vs_C1")

  seqs          <- peptable$Sequence
  uniq_mod_seqs <- unique(seqs[has_ptm])

  out_seq     <- character(0)
  out_acc     <- list()
  out_ptmpos  <- list()
  out_protptmpos <- list()
  out_ptmtype <- list()
  out_quant   <- vector("list", 0)
  out_prob   <- vector("list", 0)

  for (seq in uniq_mod_seqs) {
    mod_idx   <- which(has_ptm & seqs == seq)
    unmod_idx <- which(!has_ptm & seqs == seq)

    if (length(unmod_idx) == 0) next

    # There should be only one peptide
    if (length(unmod_idx) > 1 || length(mod_idx) > 1) {
      warning("Found more than one peptide ", seq, "!!")
      next
    }

    if ("Accession" %in% names(peptable)) {
      mod_acc <- unique_values(peptable$Accession[[mod_idx]])
      unmod_acc <- unique_values(peptable$Accession[[unmod_idx]])
      if (length(mod_acc) != 1 || length(unmod_acc) != 1 || mod_acc != unmod_acc) {
        next
      }
    }

    if ("Start" %in% names(peptable)) {
      mod_start <- unique_values(peptable$Start[[mod_idx]])
      if (length(mod_start) != 1) {
        next
      }
    }

    # Per-condition mean log2 of the same-sequence unmodified peptidoform.
    # Used only to compute occ_ref (site-specific reference occupancy).
    unmod_mean <- sapply(cond_cols, function(cols) {
      mean(as.numeric(peptable[unmod_idx, cols]), na.rm=TRUE)
    })

    # Sd across replicates averaged across conditions
    unmod_sd <- mean(sapply(cond_cols, function(cols) {
      sd(as.numeric(peptable[unmod_idx, cols]), na.rm=TRUE)
    }))


    for (mi in mod_idx) {
      # Per-condition mean log2 modified intensity
      mod_mean <- sapply(cond_cols, function(cols) {
        mean(as.numeric(peptable[mi, cols]), na.rm=TRUE)
      })
      mod_sd <- mean(sapply(cond_cols, function(cols) {
        sd(as.numeric(peptable[mi, cols]), na.rm=TRUE)
      }))

      # Protein ratio (Rprot): use only unmodified peptides from the same protein
      # accession whose sequence does NOT appear as modified (i.e., exclude any
      # counterpart unmodified peptides).  This follows requirement (c): only
      # unmodified peptides that do not have any modified version contribute to the
      # protein background ratio.  colMeans() pools rows in log2 space (geometric
      # mean in linear); mean() averages replicates within each condition.
      acc_vec        <- unlist(peptable$Accession[[mi]])
      prot_unmod_idx <- which(!has_ptm &
                                vapply(peptable$Accession,
                                       function(a) any(unlist(a) %in% acc_vec),
                                       logical(1L)) &
                                !(seqs %in% uniq_mod_seqs))
      if (length(prot_unmod_idx) == 0L) next
      prot_mean <- sapply(cond_cols, function(cols) {
        mean(colMeans(as.matrix(peptable[prot_unmod_idx, cols, drop = FALSE]), na.rm=TRUE), na.rm = TRUE)
      })
      prot_sd <- mean(sapply(cond_cols, function(cols) {
        sd(colMeans(as.matrix(peptable[prot_unmod_idx, cols, drop = FALSE]), na.rm=TRUE), na.rm = TRUE)
      }))

      log_Rm    <- mod_mean  - mod_mean[1L]   # 0 for c = 1 (reference)
      log_Ru    <- unmod_mean  - unmod_mean[1L]   # 0 for c = 1 (reference)
      log_Rprot <- prot_mean - prot_mean[1L]  # 0 for c = 1 (reference)

      Rm        <- 2^log_Rm[seq_len(NumCond-1)+1]
      Ru        <- 2^log_Ru[seq_len(NumCond-1)+1]
      Rprot     <- 2^log_Rprot[seq_len(NumCond-1)+1]

      occ    <- mean((Rprot - Ru) / (Rm - Ru))
      occ        <- c(occ, occ * Rm / Rprot)

      # Treat non-positive and NA SDs as unavailable (e.g. single replicate or
      # identical replicates produce sd = NA or sd = 0).
      if (!isTRUE(mod_sd   > 0)) mod_sd   <- NA_real_
      if (!isTRUE(unmod_sd > 0)) unmod_sd <- NA_real_
      if (!isTRUE(prot_sd  > 0)) prot_sd  <- NA_real_

      # Fill missing SDs from the available ones (×10 as uncertainty inflation).
      # Only attempt this when at least one valid SD exists; otherwise all remain
      # NA, Sigma_ab will contain NAs, and the pmvnorm block is safely skipped.
      valid_sds <- c(mod_sd, unmod_sd, prot_sd)
      valid_sds <- valid_sds[is.finite(valid_sds)]
      if (length(valid_sds) > 0) {
        ref_sd <- max(valid_sds) * 10
        if (is.na(mod_sd))   mod_sd   <- ref_sd
        if (is.na(unmod_sd)) unmod_sd <- ref_sd
        if (is.na(prot_sd))  prot_sd  <- ref_sd
      }


      # Calculating the probabilities from a multi-variate normal distirbution
      Sigma_ab <- matrix(c(
        prot_sd^2 + unmod_sd^2,   -prot_sd^2,
        -prot_sd^2,           mod_sd^2 + prot_sd^2
      ), nrow = 2, byrow = TRUE) * 2 / NumReps

      mu_ab <- cbind(log_Rprot - log_Ru, log_Rm - log_Rprot)
      p_both_neg <- p_both_pos <- vector(length=nrow(mu_ab))
      if (all(is.finite(Sigma_ab))) {
        for (i in seq_len(nrow(mu_ab))) {
          if (!any(is.na(mu_ab[i,]))) {
            p_both_neg[i] <- mvtnorm::pmvnorm(
              lower = c(-Inf, -Inf),
              upper = c(0, 0),
              mean  = mu_ab[i, ],
              sigma = Sigma_ab
            )[1]

            # P(A >= 0, B >= 0)
            p_both_pos[i] <- mvtnorm::pmvnorm(
              lower = c(0, 0),
              upper = c(Inf, Inf),
              mean  = mu_ab[i, ],
              sigma = Sigma_ab
            )[1]
          }
        }
      }

      occ_prob <- p_both_neg + p_both_pos

      occ_prob <- occ_prob[-1]

      protein_ptmpos <- NA_integer_
      if ("Start" %in% names(peptable)) {
        protein_ptmpos <- mod_start + unlist(peptable$PTMPos[[mi]]) - 1L
      }

      out_seq     <- c(out_seq, seq)
      out_acc     <- c(out_acc,     list(peptable$Accession[[mi]]))
      out_ptmpos  <- c(out_ptmpos,  list(peptable$PTMPos[[mi]]))
      out_protptmpos <- c(out_protptmpos, list(protein_ptmpos))
      out_ptmtype <- c(out_ptmtype, list(peptable$PTMType[[mi]]))
      out_quant   <- c(out_quant,   list(setNames(as.list(occ), out_cond_names)))
      out_prob    <- c(out_prob,   list(setNames(as.list(occ_prob), out_comp_names)))
    }
  }

  if (length(out_seq) == 0) {
    message("calcPTMOccupancy: no peptides with both modified and unmodified forms found.")
    return(data.frame())
  }

  quant_df               <- as.data.frame(do.call(rbind, lapply(out_quant, as.data.frame)))
  prob_df               <- as.data.frame(do.call(rbind, lapply(out_prob, as.data.frame)))
  result                 <- data.frame(Sequence = out_seq, stringsAsFactors = FALSE)
  result[out_cond_names] <- quant_df
  result[out_comp_names] <- prob_df
  result$Accession       <- out_acc
  result$PTMPos          <- out_ptmpos
  result$ProteinPTMPos   <- out_protptmpos
  result$PTMType         <- out_ptmtype

  message("  - Occupancy calculated for ", nrow(result), " modified peptidoforms across ",
          length(unique(out_seq)), " unique peptide sequences")
  result
}
