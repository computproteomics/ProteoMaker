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
#' Computes the occupancy (partial stoichiometry) for each PTM site using the
#' full three-ratio approach of Sharma et al. (2014).  For each non-reference
#' condition \eqn{c} the occupancy is:
#'
#' \deqn{occupancy_c =
#'   \frac{occ_1 \cdot R_p}{occ_1 \cdot R_p + (1 - occ_1) \cdot R_u}}
#'
#' where
#' \deqn{occ_1 = \frac{2^{\bar{l}_{p,1}}}{2^{\bar{l}_{p,1}} + 2^{\bar{l}_{u,1}}},
#'   \quad R_p = 2^{\bar{l}_{p,c} - \bar{l}_{p,1}},
#'   \quad R_u = 2^{\bar{l}_{u,c} - \bar{l}_{u,1}}}
#'
#' \eqn{occ_1} is the occupancy in the reference condition (condition 1),
#' computed within-condition from the mean log2 intensities.  \eqn{R_p} and
#' \eqn{R_u} are the between-condition intensity fold-changes for the modified
#' and unmodified peptidoforms, respectively.
#'
#' The formula can be derived directly: if \eqn{p_c} and \eqn{u_c} denote the
#' linear intensities in condition \eqn{c}, then \eqn{p_c = p_1 R_p} and
#' \eqn{u_c = u_1 R_u}, so
#' \eqn{occ_c = p_1 R_p / (p_1 R_p + u_1 R_u) = occ_1 R_p / (occ_1 R_p + (1-occ_1) R_u)}.
#' Using only \eqn{R_p} and \eqn{R_u} without \eqn{occ_1} (i.e. treating
#' \eqn{occ_1 = 0.5}) gives an incorrect result whenever the reference
#' occupancy differs from 50\%.
#'
#' Only peptide sequences observed in both a modified and an unmodified form
#' contribute to the result.  When multiple unmodified rows exist for the same
#' sequence (e.g. different charge states), their log2 intensities are pooled
#' per condition (geometric mean in linear space).  Missing values in individual
#' replicates propagate strictly: any \code{NA} makes the condition mean
#' \code{NA}, which in turn makes the occupancy \code{NA} for all conditions
#' that depend on that mean (including all non-reference conditions when the
#' \code{NA} is in the reference condition).  The reference condition column
#' (\code{C_1}) is always \code{NA} because \eqn{occ_1} is reported separately
#' as the within-condition stoichiometry, not as a ratio.
#'
#' @param peptable A data frame of peptidoform-level data as produced by the
#'   ProteoMaker pipeline (output of \code{MSRunSim} or \code{runPolySTest}).
#'   Required columns: \code{Sequence} (character), \code{PTMType} (list),
#'   \code{PTMPos} (list), \code{Accession} (list), and one column per sample
#'   as given by \code{parameters$QuantColnames}.  Quantification values must
#'   be on the log2 scale.
#' @param parameters A named list of analysis parameters.  Must contain:
#'   \describe{
#'     \item{QuantColnames}{Character vector of column names holding per-sample
#'       log2 intensities.  Samples must be ordered as all replicates of
#'       condition 1 first, then all replicates of condition 2, etc.}
#'     \item{NumCond}{Integer.  Number of experimental conditions (\eqn{\geq 2}).}
#'     \item{NumReps}{Integer.  Number of replicates per condition.}
#'   }
#'
#' @return A data frame with one row per modified peptidoform that has a
#'   matching unmodified counterpart.  Columns are:
#'   \describe{
#'     \item{Sequence}{Stripped peptide sequence.}
#'     \item{Accession}{Protein accession(s) (list column).}
#'     \item{PTMPos}{Modification positions within the peptide (list column).}
#'     \item{PTMType}{Modification types (list column).}
#'     \item{C_1, C_2, \ldots, C_NumCond}{One numeric column per condition.
#'       \code{C_1} (the reference) is always \code{NA}.  For conditions 2 and
#'       above the value is the Sharma occupancy in the range \eqn{(0, 1)}.}
#'   }
#'   Returns an empty \code{data.frame} when no modified peptides are present,
#'   when \code{NumCond < 2}, or when no modified peptide has an unmodified
#'   counterpart.
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

  # Group QuantColnames into per-condition replicate blocks
  cond_cols      <- lapply(seq_len(NumCond), function(c) {
    QuantColnames[seq.int((c - 1L) * NumReps + 1L, c * NumReps)]
  })
  out_cond_names <- paste0("C_", seq_len(NumCond))

  seqs          <- peptable$Sequence
  uniq_mod_seqs <- unique(seqs[has_ptm])

  out_seq     <- character(0)
  out_acc     <- list()
  out_ptmpos  <- list()
  out_ptmtype <- list()
  out_quant   <- vector("list", 0)

  for (seq in uniq_mod_seqs) {
    mod_idx   <- which(has_ptm & seqs == seq)
    unmod_idx <- which(!has_ptm & seqs == seq)

    if (length(unmod_idx) == 0) next

    # Per-condition mean log2 unmodified intensity.
    # colMeans() averages multiple unmodified rows in log2 space (arithmetic mean
    # in log2 = geometric mean in linear space); mean() then averages the replicates
    # within the condition.  NA in any replicate or row propagates to the condition mean.
    unmod_mean <- sapply(cond_cols, function(cols) {
      mean(colMeans(as.matrix(peptable[unmod_idx, cols, drop = FALSE])))
    })

    for (mi in mod_idx) {
      # Per-condition mean log2 modified intensity
      mod_mean <- sapply(cond_cols, function(cols) {
        mean(as.numeric(peptable[mi, cols]))
      })

      # Full Sharma et al. 3-ratio formula:
      # occ_ref from within-condition 1 (the third ratio / anchor)
      occ_ref <- 2^mod_mean[1L] / (2^mod_mean[1L] + 2^unmod_mean[1L])

      log_Rp <- mod_mean   - mod_mean[1L]    # 0 for c = 1 (reference)
      log_Ru <- unmod_mean - unmod_mean[1L]  # 0 for c = 1 (reference)
      Rp     <- 2^log_Rp
      Ru     <- 2^log_Ru
      # occ_c = (occ_ref * Rp_c) / (occ_ref * Rp_c + (1 - occ_ref) * Ru_c)
      occ        <- (occ_ref * Rp) / (occ_ref * Rp + (1 - occ_ref) * Ru)
      occ[1L]    <- NA_real_                 # reference condition is uninformative

      out_seq     <- c(out_seq, seq)
      out_acc     <- c(out_acc,     list(peptable$Accession[[mi]]))
      out_ptmpos  <- c(out_ptmpos,  list(peptable$PTMPos[[mi]]))
      out_ptmtype <- c(out_ptmtype, list(peptable$PTMType[[mi]]))
      out_quant   <- c(out_quant,   list(setNames(as.list(occ), out_cond_names)))
    }
  }

  if (length(out_seq) == 0) {
    message("calcPTMOccupancy: no peptides with both modified and unmodified forms found.")
    return(data.frame())
  }

  quant_df               <- as.data.frame(do.call(rbind, lapply(out_quant, as.data.frame)))
  result                 <- data.frame(Sequence = out_seq, stringsAsFactors = FALSE)
  result[out_cond_names] <- quant_df
  result$Accession       <- out_acc
  result$PTMPos          <- out_ptmpos
  result$PTMType         <- out_ptmtype

  message("  - Occupancy calculated for ", nrow(result), " modified peptidoforms across ",
          length(unique(out_seq)), " unique peptide sequences")
  result
}
