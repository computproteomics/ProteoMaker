#!/usr/bin/env Rscript

# Build a peptide search index from a FASTA with clear, reusable helpers.
# CLI usage:
#   Rscript inst/cmd/BuildSearchIndex.R [fasta] [enzyme] [minLen] [maxLen] [maxMC] [out.rds] [subsetN] [subsetOut]
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

# -----------------------------
# Argument parsing & parameters
# -----------------------------
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  list(
    fasta  = if (length(args) >= 1) args[[1]] else "inst/Proteomes/fasta_full_human.fasta",
    enzyme = if (length(args) >= 2) args[[2]] else "trypsin",
    minLen = if (length(args) >= 3) as.integer(args[[3]]) else 7L,
    maxLen = if (length(args) >= 4) as.integer(args[[4]]) else 30L,
    maxMC  = if (length(args) >= 5) as.integer(args[[5]]) else 2L,
    outRDS = if (length(args) >= 6) args[[6]] else NA_character_,
    subsetN = if (length(args) >= 7) as.integer(args[[7]]) else NA_integer_,
    subsetOut = if (length(args) >= 8) args[[8]] else NA_character_
  )
}

make_parameters <- function(cli) {
  list(
    PathToFasta = cli$fasta,
    Enzyme = cli$enzyme,
    PepMinLength = cli$minLen,
    PepMaxLength = cli$maxLen,
    MaxNumMissedCleavages = cli$maxMC,
    # Optional PTM defaults retained for downstream compatibility
    # PTMTypes = list(c("ph", "ox", "ac", "me", "de")),
    # ModifiableResidues = list(list(
    #   ph = c("S", "T", "Y"),
    #  ox = c("M", "W", "C", "Y"),
    #   ac = c("K"),
    #  me = c("K", "R"),
    #   de = c("D", "E")
    # ))
    PTMTypes = list(c("ox")),
    ModifiableResidues = list(list(
      # ox = c("M")
    ))
  )
}

validate_parameters <- function(p) {
  stopifnot(is.character(p$PathToFasta), nzchar(p$PathToFasta))
  if (!file.exists(p$PathToFasta)) stop("FASTA not found: ", p$PathToFasta)
  if (p$PepMinLength < 1L || p$PepMinLength > p$PepMaxLength) stop("Invalid peptide length bounds")
  if (p$MaxNumMissedCleavages < 0L) stop("maxMC must be >= 0")
  invisible(p)
}

# -----------------------------
# Core helpers (candidate to move into R/)
# -----------------------------
enzyme_regex <- function(enzyme) {
  # Mirror cleavage patterns used in R/02_Digestion.R (fastDigest)
  switch(
    enzyme,
    "trypsin"         = "(?!(RP|KP))(?=(K|R))(?!(K|R)$)",
    "trypsin.strict"  = "(?=(K|R))(?!(K|R)$)",
    "chymotrypsin.h"  = "(?!(FP|YP|PY|WP))(?=(F|Y|W))(?!(F|Y|W)$)",
    "chymotrypsin.l"  = "(?!(FP|YP|PY|WP|LP|MP))(?=(F|Y|W|L|P))(?!(F|Y|W|L|P)$)",
    "pepsin.2"        = "(?=(F|L|W|Y|A|E|Q))(?!(F|L|W|Y|A|E|Q)$)",
    "pepsin.1.3"      = "(?=(F|L))(?!(F|L)$)",
    "lysC"            = "(?=(K))(?!(K)$)",
    "argC"            = "(?!(RP))(?=(R))(?!(R)$)",
    stop(sprintf("Unsupported enzyme: %s", enzyme))
  )
}

cleavage_sites <- function(sequence, cre) {
  loc <- stringi::stri_locate_all_regex(sequence, cre)[[1]]
  if (is.null(loc) || is.na(loc[1, 1])) integer(0) else as.integer(loc[, 1])
}

segments_from_sites <- function(sequence, sites) {
  if (length(sites) == 0) {
    list(starts = 1L, stops = nchar(sequence))
  } else {
    list(starts = c(1L, sites + 1L), stops = c(sites, nchar(sequence)))
  }
}

enumerate_valid_windows <- function(starts, stops, pep_min, pep_max, max_mc) {
  valid_per_start <- vector("list", length(starts))
  n_valid <- 0L
  for (s in seq_along(starts)) {
    local <- integer(0)
    for (mc in 0:max_mc) {
      e <- s + mc
      if (e > length(stops)) break
      L <- stops[e] - starts[s] + 1L
      if (L >= pep_min && L <= pep_max) local <- c(local, mc)
    }
    valid_per_start[[s]] <- local
    n_valid <- n_valid + length(local)
  }
  list(valid_mc = valid_per_start, n_valid = n_valid)
}

normalize_IL <- function(peps) sub("I", "L", peps, perl = TRUE)

fasta_to_df <- function(path) {
  fa <- try(protr::readFASTA(file = path, legacy.mode = TRUE, seqonly = FALSE), silent = TRUE)
  if (inherits(fa, "try-error")) stop("Failed to read FASTA: ", path)
  data.frame(
    Accession = sub(".*[|]([^.]+)[|].*", "\\1", names(fa)),
    Sequence  = unlist(fa),
    stringsAsFactors = FALSE
  )
}

filter_proteins <- function(df) {
  knownAA <- c("A","L","R","K","N","M","D","F","C","P","E","S","Q","T","G","W","H","Y","I","V")
  badAA <- setdiff(LETTERS, knownAA)
  if (nrow(df) > 0) {
    keep <- !Reduce(`|`, lapply(badAA, function(x) grepl(x, df$Sequence, fixed = TRUE)))
    df <- df[keep, , drop = FALSE]
  }
  df <- df[!duplicated(df$Accession), , drop = FALSE]
  df
}

build_pep2prot_map <- function(peptides, accession, env) {
  if (length(peptides) == 0) return(invisible(NULL))
  up <- unique(peptides)
  for (p in up) {
    if (exists(p, envir = env, inherits = FALSE)) {
      cur <- get(p, envir = env, inherits = FALSE)
      if (!(accession %in% cur)) assign(p, c(cur, accession), envir = env)
    } else {
      assign(p, accession, envir = env)
    }
  }
  invisible(NULL)
}

compact_map_shared_only <- function(env) {
  if (length(ls(envir = env, all.names = TRUE)) == 0) return(list())
  keys <- ls(envir = env, all.names = TRUE)
  vals <- lapply(keys, function(k) unique(get(k, envir = env, inherits = FALSE)))
  keep <- vapply(vals, function(v) length(v) > 1, logical(1))
  if (!any(keep)) return(list())
  stats::setNames(lapply(vals[keep], function(x) sort(x)), keys[keep])
}

report_summary <- function(idx, dropped_no_windows, windows_per_protein, weights, t_secs) {
  cat(" + Proteins indexed:", length(idx), "\n")
  cat(" + Proteins with zero valid peptide windows:", dropped_no_windows, "\n")
  cat(" + Total valid peptide windows:", sum(vapply(idx, function(x) x$n_valid, integer(1))), "\n")
  cat(sprintf(" + Time: %.3f sec\n", t_secs))

  cat("\n#INDEX CONSISTENCY REPORT\n")
  cat(" + Indexed proteins:", length(idx), "\n")
  cat(" + Dropped proteins (no windows):", dropped_no_windows, "\n")
  if (length(windows_per_protein) > 0) {
    wp <- windows_per_protein[windows_per_protein > 0]
    if (length(wp) > 0) cat(" + Windows per protein (non-zero): min=", min(wp),
                            " median=", stats::median(wp),
                            " mean=", round(mean(wp), 2),
                            " max=", max(wp), "\n", sep = "")
  }
  if (length(weights) == 0) {
    cat(" + WARNING: No weights computed (empty index).\n")
  } else if (any(!is.finite(weights)) || abs(sum(weights) - 1) > 1e-6) {
    cat(" + WARNING: weights invalid (non-finite or not summing to 1).\n")
  }
}

# -----------------------------
# Peptidoform subset sampling (uniform over peptidoforms)
# -----------------------------
# Build AA -> allowed PTM types map and count per residue
build_aa_maps <- function(parameters) {
  ptm_types <- parameters$PTMTypes[[1]]
  aa_to_types <- setNames(vector("list", length = 26L), LETTERS)
  for (ptm in ptm_types) {
    aa <- parameters$ModifiableResidues[[1]][[ptm]]
    if (!is.null(aa) && length(aa) > 0) {
      for (a in aa) {
        aa_to_types[[a]] <- unique(c(aa_to_types[[a]], ptm))
      }
    }
  }
  aa_to_count <- setNames(integer(26L), LETTERS)
  for (a in names(aa_to_types)) aa_to_count[[a]] <- length(aa_to_types[[a]])
  list(aa_to_types = aa_to_types, aa_to_count = aa_to_count)
}

# Compute log number of peptidoforms for a peptide window
window_log_count <- function(seqi, start_pos, stop_pos, aa_to_count) {
  pep <- substring(seqi, start_pos, stop_pos)
  aa <- strsplit(pep, "", fixed = TRUE)[[1]]
  # choices per residue = 1 (unmodified) + number of PTM types allowed
  counts <- aa_to_count[aa]
  # use log1p to compute log(1 + x)
  sum(log1p(as.numeric(counts)))
}

# Sample a concrete peptidoform for a window uniformly over all peptidoforms
# Each residue independently picks one of (unmodified + allowed PTM types) with equal probability,
# which yields a uniform draw over the full cartesian product of choices.
sample_uniform_peptidoform_in_window <- function(seqi, start_pos, stop_pos, aa_to_types) {
  pep <- substring(seqi, start_pos, stop_pos)
  aa <- strsplit(pep, "", fixed = TRUE)[[1]]
  ptm_positions <- integer(0)
  ptm_types <- character(0)
  for (j in seq_along(aa)) {
    allowed <- aa_to_types[[aa[j]]]
    k <- length(allowed)
    # draw among 0..k uniformly; 0=no PTM here
    draw <- sample.int(k + 1L, 1L) - 1L
    if (draw > 0L) {
      chosen <- allowed[[draw]]
      ptm_positions <- c(ptm_positions, j)
      ptm_types <- c(ptm_types, chosen)
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

# Sample N peptidoforms uniformly from the full space without enumerating all combinations materialized
sample_peptidoforms_from_index <- function(index, parameters, N) {
  if (is.na(N) || N <= 0) return(NULL)
  maps <- build_aa_maps(parameters)
  # Build window table with log counts
  pidxs <- integer(0); sidxs <- integer(0); mcs <- integer(0)
  starts <- integer(0); stops <- integer(0); logc <- numeric(0)
  for (pi in seq_along(index$proteins)) {
    e <- index$proteins[[pi]]
    if (is.null(e)) next
    if (length(e$valid_mc) == 0) next
    s_rep <- rep(seq_along(e$starts), lengths(e$valid_mc))
    MC <- unlist(e$valid_mc)
    if (length(MC) == 0) next
    st_pos <- e$starts[s_rep]
    en_pos <- e$stops[s_rep + MC]
    # compute log counts for each peptide efficiently
    lc <- mapply(function(a,b) window_log_count(e$sequence, a, b, maps$aa_to_count), st_pos, en_pos)
    pidxs <- c(pidxs, rep.int(pi, length(MC)))
    sidxs <- c(sidxs, s_rep)
    mcs   <- c(mcs, MC)
    starts<- c(starts, st_pos)
    stops <- c(stops, en_pos)
    logc  <- c(logc, lc)
  }
  if (length(logc) == 0) return(NULL)
  # stable probabilities from log counts
  maxl <- max(logc)
  probs <- exp(logc - maxl)
  probs <- probs / sum(probs)
  # sample indices with replacement to allow N > unique count
  pick <- sample.int(length(probs), size = N, replace = TRUE, prob = probs)
  # materialize chosen peptidoforms
  out <- vector("list", length(pick))
  for (i in seq_along(pick)) {
    k <- pick[[i]]
    e <- index$proteins[[pidxs[[k]]]]
    sp <- starts[[k]]; ep <- stops[[k]]
    spf <- sample_uniform_peptidoform_in_window(e$sequence, sp, ep, maps$aa_to_types)
    out[[i]] <- data.frame(
      Accession = e$accession,
      Peptide = spf$Peptide,
      Start = spf$Start,
      Stop = spf$Stop,
      MC = mcs[[k]],
      stringsAsFactors = FALSE
    )
    out[[i]]$PTMPos <- spf$PTMPos
    out[[i]]$PTMType <- spf$PTMType
  }
  do.call(rbind, out)
}

# Render a human-readable peptidoform by annotating PTMs in the sequence
# Example: S at position 3 with type 'ph' becomes 'S[ph]'
annotate_peptidoform <- function(pep, pos, typ) {
  if (length(pos) == 0 || length(unlist(pos)) == 0) return(pep)
  p <- unlist(pos)
  t <- unlist(typ)
  chars <- strsplit(pep, "", fixed = TRUE)[[1]]
  # Insert annotations after each modified residue (process from right to left to keep indices stable)
  ord <- order(p, decreasing = TRUE)
  for (i in ord) {
    k <- p[[i]]
    if (k >= 1 && k <= length(chars)) {
      chars[k] <- paste0(chars[k], "[", t[[i]], "]")
    }
  }
  paste0(chars, collapse = "")
}

# -----------------------------
# Index builder
# -----------------------------
build_search_index <- function(parameters) {
  cat("\n#INDEX BUILD - Start\n")
  cat(" + FASTA:", parameters$PathToFasta, "\n")
  cat(" + Enzyme:", parameters$Enzyme,
      " minLen:", parameters$PepMinLength,
      " maxLen:", parameters$PepMaxLength,
      " maxMC:", parameters$MaxNumMissedCleavages, "\n")

  t0 <- proc.time()[[3]]
  df <- filter_proteins(fasta_to_df(parameters$PathToFasta))
  cat(" + Proteins after filtering:", nrow(df), "\n")
  if (nrow(df) == 0) stop("No proteins left after filtering")

  cre <- enzyme_regex(parameters$Enzyme)

  idx <- vector("list", nrow(df))
  dropped_no_windows <- 0L
  windows_per_protein <- integer(nrow(df))
  pep2prot_env <- new.env(parent = emptyenv())

  for (i in seq_len(nrow(df))) {
    seqi <- df$Sequence[i]
    sites <- cleavage_sites(seqi, cre)
    seg   <- segments_from_sites(seqi, sites)
    win   <- enumerate_valid_windows(seg$starts, seg$stops,
                                     parameters$PepMinLength,
                                     parameters$PepMaxLength,
                                     parameters$MaxNumMissedCleavages)

    windows_per_protein[i] <- win$n_valid
    if (win$n_valid == 0L) {
      idx[[i]] <- NULL
      dropped_no_windows <- dropped_no_windows + 1L
      next
    }

    # Optional PTM positions (kept for downstream consumers)
    modpos <- list()
    if (!is.null(parameters$ModifiableResidues) && length(parameters$ModifiableResidues) > 0) {
      ptms <- parameters$PTMTypes[[1]]
      for (ptm in ptms) {
        aa <- parameters$ModifiableResidues[[1]][[ptm]]
        modpos[[ptm]] <- if (!is.null(aa)) which(strsplit(seqi, "")[[1]] %in% aa) else integer(0)
      }
    }

    idx[[i]] <- list(
      accession = df$Accession[i],
      sequence  = seqi,
      starts    = seg$starts,
      stops     = seg$stops,
      valid_mc  = win$valid_mc,
      n_valid   = win$n_valid,
      modpos    = modpos
    )

    # Build peptide -> protein map (I->L normalized)
    starts_idx <- rep(seq_along(seg$starts), lengths(win$valid_mc))
    MC <- unlist(win$valid_mc)
    if (length(MC) > 0) {
      stops_idx <- starts_idx + MC
      pep <- substring(seqi, seg$starts[starts_idx], seg$stops[stops_idx])
      build_pep2prot_map(normalize_IL(pep), df$Accession[i], pep2prot_env)
    }
  }

  # Compact index and compute weights
  keep_mask <- !vapply(idx, is.null, logical(1))
  idx <- idx[keep_mask]
  raw_w <- if (length(idx) > 0) vapply(idx, function(x) x$n_valid, integer(1)) else integer(0)
  weights <- if (length(raw_w) == 0) {
    numeric(0)
  } else if (sum(raw_w) > 0) {
    raw_w / sum(raw_w)
  } else {
    rep(1 / length(raw_w), length(raw_w))
  }

  t1 <- proc.time()[[3]]
  report_summary(idx, dropped_no_windows, windows_per_protein, weights, t1 - t0)

  pep2prot <- compact_map_shared_only(pep2prot_env)
  if (length(pep2prot) > 0) {
    n_shared <- length(pep2prot)
    cat(" + Shared peptides (pep2prot entries):", n_shared, "\n")
    amb_deg <- vapply(pep2prot, length, integer(1))
    cat(" + Peptide sharedness (proteins per peptide): min=", min(amb_deg),
        " median=", stats::median(amb_deg),
        " mean=", round(mean(amb_deg), 2),
        " max=", max(amb_deg), "\n", sep = "")
    top_idx <- order(amb_deg, decreasing = TRUE)[seq_len(min(5, length(amb_deg)))]
    cat(" + Top shared peptides (truncated):\n")
    for (k in top_idx) {
      nm <- names(pep2prot)[k]
      cat("    - ", substr(nm, 1, 30), if (nchar(nm) > 30) "..." else "",
          " -> ", amb_deg[k], " proteins\n", sep = "")
    }
  }

  res <- list(
    proteins = idx,
    weights = weights,
    pep2prot = pep2prot,
    params = parameters,
    build_time_sec = t1 - t0
  )
}

# -----------------------------
# Entrypoint
# -----------------------------
cli <- parse_args()
params <- validate_parameters(make_parameters(cli))
res <- build_search_index(params)

if (!is.na(cli$outRDS)) {
  saveRDS(res, file = cli$outRDS)
  cat("\nSaved index to:", cli$outRDS, "\n")
}

# Optional: sample a random set of peptidoforms uniformly across the full space
if (!is.na(cli$subsetN) && cli$subsetN > 0) {
  cat("\n#PEPTIDOFORMS SUBSET - Start\n")
  subset <- sample_peptidoforms_from_index(res, params, cli$subsetN)
  if (is.null(subset)) {
    cat(" + No peptidoforms could be sampled.\n")
  } else {
    cat(" + Sampled peptidoforms:", nrow(subset), "(uniform across peptidoforms)\n")
    headn <- min(5L, nrow(subset))
    cat(" + Preview (", headn, ") [PTMs in brackets]:\n", sep = "")
    prev <- utils::head(subset, headn)
    prev$Peptidoform <- mapply(annotate_peptidoform, prev$Peptide, prev$PTMPos, prev$PTMType, SIMPLIFY = TRUE)
    print(prev[, c("Accession", "Peptidoform", "Start", "Stop", "MC")], row.names = FALSE)
    if (!is.na(cli$subsetOut) && nzchar(cli$subsetOut)) {
      saveRDS(subset, file = cli$subsetOut)
      cat(" + Saved peptidoform subset to:", cli$subsetOut, "\n")
    }
  }
  cat("#PEPTIDOFORMS SUBSET - Finish\n")
}
