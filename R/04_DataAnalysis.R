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
#' Computes the occupancy (partial stoichiometry) for each PTM site from
#' peptidoform-level quantification data. Occupancy is defined as the fraction
#' of protein molecules that carry the modification at a given site:
#'
#' \deqn{occupancy = \frac{I_{mod}}{I_{mod} + I_{unmod}}}
#'
#' where \eqn{I_{mod}} and \eqn{I_{unmod}} are the linear-scale intensities of
#' the modified and unmodified peptidoforms, respectively.  The input
#' intensities are assumed to be on the log2 scale, as produced by the
#' ProteoMaker pipeline.
#'
#' Only peptide sequences that are observed in both a modified and an
#' unmodified form contribute to the result.  When multiple unmodified rows
#' exist for the same sequence (e.g., from different charge states), their
#' linear-scale intensities are averaged before computing occupancy.
#' Missing values (\code{NA}) in individual samples are propagated: if either
#' the modified or unmodified intensity is \code{NA} for a sample, the
#' occupancy for that sample is also \code{NA}.  When multiple unmodified rows
#' are present and any one of them has \code{NA} for a given sample, the
#' averaged unmodified signal for that sample is also \code{NA}, so occupancy
#' is \code{NA} as well.
#'
#' @param peptable A data frame of peptidoform-level data as produced by the
#'   ProteoMaker pipeline (output of \code{MSRunSim} or \code{runPolySTest}).
#'   Required columns: \code{Sequence} (character), \code{PTMType} (list),
#'   \code{PTMPos} (list), \code{Accession} (list), and one column per sample
#'   as given by \code{parameters$QuantColnames}.  Quantification values must
#'   be on the log2 scale.
#' @param parameters A named list of analysis parameters.  Must contain at
#'   least \code{QuantColnames}, a character vector of the column names that
#'   hold the per-sample log2 intensities.
#'
#' @return A data frame with one row per modified peptidoform that has a
#'   matching unmodified counterpart.  Columns are:
#'   \describe{
#'     \item{Sequence}{Stripped peptide sequence.}
#'     \item{Accession}{Protein accession(s) (list column).}
#'     \item{PTMPos}{Modification positions within the peptide (list column).}
#'     \item{PTMType}{Modification types (list column).}
#'     \item{<sample>}{One numeric column per entry in \code{QuantColnames},
#'       containing the per-sample occupancy in the range \eqn{[0, 1]}.}
#'   }
#'   Returns an empty \code{data.frame} when no modified peptides are present
#'   or when no modified peptide has an unmodified counterpart.
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

  has_ptm <- lengths(peptable$PTMType) > 0

  if (sum(has_ptm) == 0) {
    message("calcPTMOccupancy: no modified peptides found; returning empty table.")
    return(data.frame())
  }

  message(" + Calculating PTM site occupancy")

  seqs <- peptable$Sequence
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

    # Unmodified linear-scale signal; average across multiple unmodified rows
    unmod_mat <- 2^as.matrix(peptable[unmod_idx, QuantColnames, drop = FALSE])
    unmod_lin <- colMeans(unmod_mat)

    for (mi in mod_idx) {
      mod_lin <- 2^as.numeric(peptable[mi, QuantColnames])
      total   <- mod_lin + unmod_lin
      occ     <- ifelse(is.na(mod_lin) | is.na(unmod_lin) | total == 0, NA_real_, mod_lin / total)

      out_seq     <- c(out_seq, seq)
      out_acc     <- c(out_acc,     list(peptable$Accession[[mi]]))
      out_ptmpos  <- c(out_ptmpos,  list(peptable$PTMPos[[mi]]))
      out_ptmtype <- c(out_ptmtype, list(peptable$PTMType[[mi]]))
      out_quant   <- c(out_quant,   list(setNames(as.list(occ), QuantColnames)))
    }
  }

  if (length(out_seq) == 0) {
    message("calcPTMOccupancy: no peptides with both modified and unmodified forms found.")
    return(data.frame())
  }

  quant_df           <- as.data.frame(do.call(rbind, lapply(out_quant, as.data.frame)))
  result             <- data.frame(Sequence = out_seq, stringsAsFactors = FALSE)
  result[QuantColnames] <- quant_df
  result$Accession   <- out_acc
  result$PTMPos      <- out_ptmpos
  result$PTMType     <- out_ptmtype

  message("  - Occupancy calculated for ", nrow(result), " modified peptidoforms across ",
          length(unique(out_seq)), " unique peptide sequences")
  result
}
