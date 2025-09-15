#!/usr/bin/env Rscript

# Standalone script to build and time a peptide search index from a FASTA.
# Usage:
#   Rscript inst/cmd/BuildSearchIndex.R [fasta] [enzyme] [minLen] [maxLen] [maxMC] [out.rds]
# Defaults:
#   fasta   = inst/Proteomes/fasta_full_human.fasta
#   enzyme  = trypsin
#   minLen  = 7
#   maxLen  = 30
#   maxMC   = 2
#   out.rds = (optional) path to save the index via saveRDS

suppressPackageStartupMessages({
  requireNamespace("protr", quietly = TRUE)
  requireNamespace("stringi", quietly = TRUE)
})

args <- commandArgs(trailingOnly = TRUE)
fasta <- if (length(args) >= 1) args[[1]] else "inst/Proteomes/fasta_full_human.fasta"
enzyme <- if (length(args) >= 2) args[[2]] else "trypsin"
minLen <- if (length(args) >= 3) as.integer(args[[3]]) else 7L
maxLen <- if (length(args) >= 4) as.integer(args[[4]]) else 30L
maxMC  <- if (length(args) >= 5) as.integer(args[[5]]) else 2L
outRDS <- if (length(args) >= 6) args[[6]] else NA_character_

# Minimal parameter list needed for indexing
parameters <- list(
  PathToFasta = fasta,
  Enzyme = enzyme,
  PepMinLength = minLen,
  PepMaxLength = maxLen,
  MaxNumMissedCleavages = maxMC,
  # Provide defaults for PTM definitions (not used for index speed, but kept for completeness)
  PTMTypes = list(c("ph", "ox", "ac", "me", "de")),
  ModifiableResidues = list(list(
    ph = c("S", "T", "Y"),
    ox = c("M", "W", "C", "Y"),
    ac = c("K"),
    me = c("K", "R"),
    de = c("D", "E")
  ))
)

enzyme_regex <- function(enzyme) {
  # Mirror the cleavage patterns used in R/02_Digestion.R (fastDigest)
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

build_search_index <- function(parameters) {
  cat("\n#INDEX BUILD - Start\n")
  cat(" + FASTA:", parameters$PathToFasta, "\n")
  cat(" + Enzyme:", parameters$Enzyme, " minLen:", parameters$PepMinLength,
      " maxLen:", parameters$PepMaxLength, " maxMC:", parameters$MaxNumMissedCleavages, "\n")

  t0 <- proc.time()[[3]]
  fasta_obj <- try(protr::readFASTA(file = parameters$PathToFasta, legacy.mode = TRUE, seqonly = FALSE), silent = TRUE)
  if (inherits(fasta_obj, "try-error")) stop("Failed to read FASTA: ", parameters$PathToFasta)
  df <- data.frame(
    Accession = sub(".*[|]([^.]+)[|].*", "\\1", names(fasta_obj)),
    Sequence  = unlist(fasta_obj),
    stringsAsFactors = FALSE
  )
  # Filter unusual amino acids and duplicates
  knownAA <- c("A","L","R","K","N","M","D","F","C","P","E","S","Q","T","G","W","H","Y","I","V")
  badAA <- setdiff(LETTERS, knownAA)
  if (nrow(df) > 0) {
    keep <- !Reduce(`|`, lapply(badAA, function(x) grepl(x, df$Sequence)))
    df <- df[keep, , drop = FALSE]
  }
  df <- df[!duplicated(df$Accession), , drop = FALSE]

  cat(" + Proteins after filtering:", nrow(df), "\n")
  if (nrow(df) == 0) stop("No proteins left after filtering")

  cre <- enzyme_regex(parameters$Enzyme)

  # Precompute per-protein cleavage segments and valid windows per start
  idx <- vector("list", length = nrow(df))
  total_valid <- 0L
  # Peptide (normalized I->L) -> proteins mapping; we fill only for peptides observed in >= 2 proteins later
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
    cs <- stringi::stri_locate_all_regex(seq, cre)[[1]][,1]
    if (length(cs) == 0 || is.na(cs[1])) {
      # No cleavage sites -> single peptide candidate if length fits
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
      total_valid <- total_valid + n_valid
      # Optional: precompute modifiable positions per PTM type (can be used downstream)
      modpos <- list()
      if (!is.null(parameters$ModifiableResidues) && length(parameters$ModifiableResidues) > 0) {
        ptms <- parameters$PTMTypes[[1]]
        for (ptm in ptms) {
          aa <- parameters$ModifiableResidues[[1]][[ptm]]
          if (!is.null(aa)) {
            modpos[[ptm]] <- which(strsplit(seq, "")[[1]] %in% aa)
          } else {
            modpos[[ptm]] <- integer(0)
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
  cat(" + Proteins indexed:", length(idx), "\n")
  cat(" + Total valid peptide windows:", sum(vapply(idx, function(x) x$n_valid, integer(1))), "\n")
  cat(sprintf(" + Time: %.3f sec\n", t1 - t0))
  cat("#INDEX BUILD - Finish\n")

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

res <- build_search_index(parameters)

if (!is.na(outRDS)) {
  saveRDS(res, file = outRDS)
  cat("\nSaved index to:", outRDS, "\n")
}
