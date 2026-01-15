################################################################################
#               Data analysis of the outputs of the peptide quan.              #
################################################################################

#' Summarize Protein Abundance from Peptide-Level Data
#'
#' This function summarizes protein abundance from peptide-level data. It removes modified peptides, groups peptides by their parent proteins, and applies a summarization method (e.g., sum of top 3, median polish) to estimate protein abundance. The function supports parallel processing to speed up the summarization process.
#'
#' @param peptable A data frame containing peptide-level data, including accession numbers, sequence, PTM types, and quantification columns.
#' @param parameters A list of parameters, including:
#' \describe{
#'   \item{ProtSummarization}{The method used for protein summarization, e.g., "sum.top3" or "medpolish".}
#'   \item{MinUniquePep}{The minimum number of unique peptides required to summarize a protein.}
#'   \item{QuantColnames}{The names of the columns containing quantification data.}
#'   \item{Cores}{The number of cores to use for parallel processing.}
#'   \item{ClusterType}{The type of cluster to use for parallel processing (e.g., "FORK", "PSOCK").}
#' }
#'
#' @return A data frame containing summarized protein-level data, where each row represents a protein, and the columns include protein information and summarized quantification data.
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
    # Remove all modified peptides
    peptable <- peptable[!sapply(peptable$PTMType, function(x) length(x) > 0), ]
    message("  - Removed all modified peptides, remaining number of peptides: ", nrow(peptable), "")
  } else {
    message("  - Keeping modified peptides")
  }
  message("  - Remaining number of peptides: ", nrow(peptable), "")

  message(" + Protein summarisation")

  message("  - Number of minimum unique peptide per protein: ", minUniquePep, "")
  message("  - Protein summarisation using the ", method, " approach.")

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
  message("  - Summarizing proteins, this can take a while")

  # Function to summarize protein groups
  summarizeProtein <- function(tmp) {
    out <- NULL
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
          summed <- tmp - median(summed, na.rm = TRUE)
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
          fit <- rlm(intensity ~ peptide + sample, data = long_df, na.action = na.omit)
          coefs <- coef(fit)
          sample_coefs <- coefs[grep("^sample", names(coefs))]

          # Include intercept (baseline)
          baseline <- coefs["(Intercept)"]
          sample_names <- gsub("^sample", "", names(sample_coefs))

          # Set full abundance vector
          out <- setNames(baseline + sample_coefs, sample_names)
        } else {
          out <- tmp[1, QuantColnames]
        }
      } else {
        stop("No valid method provided!")
      }
    }
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
      out <- tmp[1, ]
      tout <- summarizeProtein(tmp[, QuantColnames, drop = F])
      if (!is.null(tout)) {
        out[QuantColnames] <- tout
        # add other information
        out[other_cols] <- sapply(tmp[, other_cols], function(x) paste(unlist(x), collapse = ";"))
        # ensure unique row name per protein group (use group key)
        rownames(out) <- names(prot_ind)[i]
      } else {
        out <- NULL
      }
      return(out)
    })
    parallel::stopCluster(cluster)
  } else {
    proteins <- lapply(1:(length(prot_ind) - 1), function(i) {
      tmp <- as.data.frame(peptable[prot_ind[i]:(prot_ind[i + 1] - 1), ])
      key <- if ("Peptidoform" %in% names(tmp)) tmp$Peptidoform else tmp$Sequence
      rownames(tmp) <- make.unique(key)
      out <- tmp[1, ]
      tout <- summarizeProtein(tmp[, QuantColnames, drop = F])
      if (!is.null(tout)) {
        out[QuantColnames] <- tout
        # add other information
        out[other_cols] <- sapply(tmp[, other_cols], function(x) paste(unlist(x), collapse = ";"))
        rownames(out) <- names(prot_ind)[i]
      } else {
        out <- NULL
      }
      return(out)
    })
  }
  # join all protein data
  proteins <- Filter(Negate(is.null), proteins)
  protmat <- do.call(rbind, proteins)
  protmat[protmat == -Inf] <- NA
  protmat <- protmat[rowSums(is.na(protmat[, QuantColnames])) < length(QuantColnames), ]

  cat("  - Finished summarizing into ", nrow(protmat), " proteins")
  #  for (i in parameters$QuantColnames) protmat[,i] <- as.numeric(protmat[,i])

  return(protmat)
}
