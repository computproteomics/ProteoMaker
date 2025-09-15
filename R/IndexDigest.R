################################################################################
#                         INDEX-ACCELERATED DIGESTION                           #
################################################################################

#' Build peptide template from indexed protein
#'
#' Given a single protein entry from the SearchIndex (as built by
#' inst/cmd/BuildSearchIndex.R or a compatible builder), expand all valid
#' (start, mc) windows into a peptide table matching fastDigest output
#' (without PTMs) and including precomputed M/Z values. This template can be
#' reused across all proteoforms of the same protein to avoid repeated regex
#' cleavage and mass calculations.
#'
#' @param idx_protein A single element from parameters$SearchIndex$proteins.
#' @return A data.frame with columns Peptide, Start, Stop, MC, MZ1, MZ2, MZ3.
#' @keywords internal
protein_peptide_template_from_index <- function(idx_protein) {

  AA.mass <- c(
    A = 71.03711, R = 156.10111, N = 114.04293, D = 115.02694, C = 103.00919,
    E = 129.04259, Q = 128.05858, G = 57.02146, H = 137.05891, I = 113.08406,
    L = 113.08406, K = 128.09496, M = 131.04049, F = 147.06841, P = 97.05276,
    S = 87.03203, T = 101.04768, W = 186.07931, Y = 163.06333, V = 99.06841
  )

  seq <- idx_protein$sequence

  # Expand valid windows
  starts_idx <- rep(seq_along(idx_protein$starts), lengths(idx_protein$valid_mc))
  MC <- unlist(idx_protein$valid_mc)
  stops_idx <- starts_idx + MC
  start_pos <- idx_protein$starts[starts_idx]
  stop_pos  <- idx_protein$stops[stops_idx]

  pep <- substring(seq, start_pos, stop_pos)

  # Compute monoisotopic peptide mass + water, then M/Z for +1, +2, +3
  pep_mass <- vapply(strsplit(pep, ""), function(aas) sum(AA.mass[aas]) + 18.01528, numeric(1))
  MZ1 <- pep_mass + 1.007276466
  MZ2 <- (pep_mass + 2 * 1.007276466) / 2
  MZ3 <- (pep_mass + 3 * 1.007276466) / 3

  data.frame(
    Peptide = pep,
    Start = start_pos,
    Stop = stop_pos,
    MC = MC,
    MZ1 = MZ1, MZ2 = MZ2, MZ3 = MZ3,
    stringsAsFactors = FALSE
  )
}

#' Find matching index entry for a proteoform
#'
#' Locate the corresponding indexed protein entry for a given proteoform row by
#' matching Accession and Sequence.
#'
#' @param proteoform Single-row data.frame with Accession and Sequence fields.
#' @param index List as in parameters$SearchIndex.
#' @return The matching index entry or NULL if not found.
#' @keywords internal
find_index_for_proteoform <- function(proteoform, index) {
  if (is.null(index) || is.null(index$proteins)) return(NULL)
  for (j in seq_along(index$proteins)) {
    ip <- index$proteins[[j]]
    if (identical(ip$accession, proteoform$Accession) && identical(ip$sequence, proteoform$Sequence)) {
      return(ip)
    }
  }
  return(NULL)
}

#' Build search index for fast digestion and wrong-ID sampling
#'
#' Precompute per-protein cleavage segments and valid peptide windows given
#' enzyme rules, peptide length bounds, and max missed cleavages. Also records
#' positions of modifiable residues per PTM type. The index can be attached to
#' the parameters list as `SearchIndex` and used by `proteoformDigestion`.
#'
#' @param parameters List containing at least PathToFasta, Enzyme, PepMinLength,
#' PepMaxLength, MaxNumMissedCleavages, PTMTypes, ModifiableResidues.
#' @return A list with fields: proteins (list of entries), weights (per protein),
#' params (subset of parameters), build_time_sec.
#' @keywords internal
## build_search_index is defined in 01_GenerateGroundTruth.R to keep ground-truth
## related functionality together.

#' Build a digestion search index from an in-memory protein data frame
#'
#' This helper constructs the same index structure as the standalone
#' builder but operates on a data.frame with columns Accession and Sequence.
#' It mirrors the enzyme cleavage rules used by fastDigest to ensure
#' consistency in peptide windows.
#'
#' @param df Data frame with columns Accession and Sequence (one row per protein).
#' @param parameters List containing Enzyme, PepMinLength, PepMaxLength, MaxNumMissedCleavages,
#' and optionally PTMTypes and ModifiableResidues to annotate modifiable positions.
#' @return A list with fields: proteins (list), weights (numeric vector), params (list), build_time_sec (numeric).
#' @keywords internal
build_search_index_from_df <- function(df, parameters) {
  t0 <- proc.time()[[3]]

  enzyme_regex <- function(enzyme) {
    if (enzyme == "trypsin") return("(?!(RP|KP))(?=(K|R))(?!(K|R)$)")
    if (enzyme == "trypsin.strict") return("(?=(K|R))(?!(K|R)$)")
    if (enzyme == "chymotrypsin.h") return("(?!(FP|YP|PY|WP))(?=(F|Y|W))(?!(F|Y|W)$)")
    if (enzyme == "chymotrypsin.l") return("(?!(FP|YP|PY|WP|LP|MP))(?=(F|Y|W|L|P))(?!(F|Y|W|L|P)$)")
    if (enzyme == "pepsin.2") return("(?=(F|L|W|Y|A|E|Q))(?!(F|L|W|Y|A|E|Q)$)")
    if (enzyme == "pepsin.1.3") return("(?=(F|L))(?!(F|L)$)")
    if (enzyme == "lysC") return("(?=(K))(?!(K)$)")
    if (enzyme == "argC") return("(?!(RP))(?=(R))(?!(R)$)")
    stop(sprintf("Unsupported enzyme: %s", enzyme))
  }

  cre <- enzyme_regex(parameters$Enzyme)

  idx <- vector("list", length = nrow(df))
  # Peptide (normalized I->L) -> proteins mapping; collect and compact to shared only at the end
  pep2prot_env <- new.env(parent = emptyenv())
  add_map <- function(pep_vec, acc) {
    if (length(pep_vec) == 0) return(invisible(NULL))
    up <- unique(pep_vec)
    for (p in up) {
      if (exists(p, envir = pep2prot_env, inherits = FALSE)) {
        current <- get(p, envir = pep2prot_env, inherits = FALSE)
        if (!(acc %in% current)) assign(p, c(current, acc), envir = pep2prot_env)
      } else {
        assign(p, acc, envir = pep2prot_env)
      }
    }
  }

  for (i in seq_len(nrow(df))) {
    seq <- df$Sequence[i]
    cs <- stringi::stri_locate_all_regex(seq, cre)[[1]][, 1]
    if (length(cs) == 0 || is.na(cs[1])) {
      starts <- 1L
      stops  <- nchar(seq)
    } else {
      starts <- c(1L, cs + 1L)
      stops  <- c(cs, nchar(seq))
    }

    valid_per_start <- vector("list", length(starts))
    n_valid <- 0L
    for (s in seq_along(starts)) {
      local <- integer(0)
      for (mc in 0:parameters$MaxNumMissedCleavages) {
        e <- s + mc
        if (e > length(stops)) break
        L <- stops[e] - starts[s] + 1L
        if (L >= parameters$PepMinLength && L <= parameters$PepMaxLength) local <- c(local, mc)
      }
      valid_per_start[[s]] <- local
      n_valid <- n_valid + length(local)
    }

    if (n_valid == 0L) {
      idx[[i]] <- NULL
    } else {
      # Optional: modifiable positions per PTM type
      modpos <- list()
      if (!is.null(parameters$ModifiableResidues) && length(parameters$ModifiableResidues) > 0) {
        ptms <- parameters$PTMTypes[[1]]
        if (!is.null(ptms)) {
          for (ptm in ptms) {
            aa <- parameters$ModifiableResidues[[1]][[ptm]]
            if (!is.null(aa)) {
              modpos[[ptm]] <- which(strsplit(seq, "")[[1]] %in% aa)
            } else {
              modpos[[ptm]] <- integer(0)
            }
          }
        }
      }
      idx[[i]] <- list(
        accession = df$Accession[i],
        sequence  = seq,
        starts    = starts,
        stops     = stops,
        valid_mc  = valid_per_start,
        n_valid   = n_valid,
        modpos    = modpos
      )

      # Update peptide->protein mapping for this protein
      starts_idx <- rep(seq_along(starts), lengths(valid_per_start))
      MC <- unlist(valid_per_start)
      if (length(MC) > 0) {
        stops_idx <- starts_idx + MC
        start_pos <- starts[starts_idx]
        stop_pos  <- stops[stops_idx]
        pep <- substring(seq, start_pos, stop_pos)
        pep_norm <- gsub("[I]", "L", pep)
        add_map(pep_norm, df$Accession[i])
      }
    }
  }

  idx <- idx[!vapply(idx, is.null, logical(1))]
  weights <- vapply(idx, function(x) x$n_valid, integer(1))
  weights <- weights / sum(weights)

  t1 <- proc.time()[[3]]
  # Compact peptide->protein mapping: keep only peptides shared by >= 2 proteins
  if (length(ls(envir = pep2prot_env, all.names = TRUE)) > 0) {
    keys <- ls(envir = pep2prot_env, all.names = TRUE)
    vals <- lapply(keys, function(k) unique(get(k, envir = pep2prot_env, inherits = FALSE)))
    keep <- vapply(vals, function(v) length(v) > 1, logical(1))
    if (any(keep)) {
      pep2prot <- setNames(lapply(vals[keep], function(x) sort(x)), keys[keep])
    } else {
      pep2prot <- list()
    }
  } else {
    pep2prot <- list()
  }

  list(proteins = idx, weights = weights, pep2prot = pep2prot, params = parameters, build_time_sec = t1 - t0)
}

#' Build a simple protein index from Accession/Sequence (no enzyme settings required)
#'
#' This index only stores per-protein accession and sequence plus simple weights.
#' It can be used at MSRun time to generate decoy peptide templates on-the-fly
#' using fastDigest and the current digestion parameters.
#'
#' @param df Data frame with columns Accession and Sequence.
#' @return A list with fields: proteins (list of list(accession, sequence)), weights (uniform), params = NULL.
#' @keywords internal
build_protein_index_from_df <- function(df) {
  idx <- lapply(seq_len(nrow(df)), function(i) list(accession = df$Accession[i], sequence = df$Sequence[i]))
  weights <- rep(1/length(idx), length(idx))
  list(proteins = idx, weights = weights, params = NULL, build_time_sec = NA_real_)
}
