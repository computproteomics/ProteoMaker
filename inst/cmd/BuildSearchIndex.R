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
    PTMTypes = list(c("ph", "ox", "ac", "me", "de")),
    ModifiableResidues = list(list(
      ph = c("S", "T", "Y"),
     ox = c("M", "W", "C", "Y"),
      ac = c("K"),
     me = c("K", "R"),
      de = c("D", "E")
    ))
    # PTMTypes = list(c("ph","ox")),
    # ModifiableResidues = list(list(
    #     ph = c("S", "T", "Y"),
    #   ox = c("M")
    # ))
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

normalize_IL <- function(peps) gsub("I", "L", peps, perl = TRUE)

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
  if (!is.character(peptides)) peptides <- as.character(peptides)
  bad <- which(is.na(peptides) | peptides == "")
  if (length(bad) > 0) {
    cat(" + WARNING: encountered ", length(bad), " invalid peptide names (empty/NA) for accession ", accession, "; showing up to 5 examples:\n", sep = "")
    print(utils::head(peptides[bad], 5))
    peptides <- peptides[-bad]
    if (length(peptides) == 0) return(invisible(NULL))
  }
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

sample_peptidoforms <- function(index, parameters, N) {
  if (is.null(index) || is.null(index$proteins) || N <= 0) return(NULL)
  # Build PTM assignment map once
  maps <- build_aa_maps(parameters)

  # Protein weights: pf_total (sum of per-window counts). If all zero (shouldn't), fall back to n_valid
  pf <- vapply(index$proteins, function(e) if (is.null(e)) 0 else as.numeric(e$pf_total), numeric(1))
  if (!any(pf > 0)) pf <- vapply(index$proteins, function(e) if (is.null(e)) 0L else e$n_valid, integer(1))
  total_pf <- sum(pf)
  if (total_pf <= 0) {
    message("[subset] No peptidoforms available: total_pf=0; check length/missed-cleavage constraints.")
    return(NULL)
  }

  # Draw proteins for each peptidoform
  prot_indices <- sample.int(length(index$proteins), size = N, replace = TRUE, prob = pf)
  # Group draws per protein to amortize per-protein work
  by_prot <- split(seq_along(prot_indices), prot_indices)

  out <- vector("list", N)
  out_ptr <- 1L
  for (pi_chr in names(by_prot)) {
    pi <- as.integer(pi_chr)
    e <- index$proteins[[pi]]
    idxs <- by_prot[[pi_chr]]
    n_draws <- length(idxs)
    # Derive per-window weights within this protein
    if (!is.null(e$win_count) && length(e$win_count) > 0 && sum(e$win_count) > 0 && all(is.finite(e$win_count))) {
      win_probs <- e$win_count / sum(e$win_count)
      pick_win <- sample.int(length(win_probs), size = n_draws, replace = TRUE, prob = win_probs)
      for (jj in seq_len(n_draws)) {
        k <- pick_win[[jj]]
        sp <- e$win_start[[k]]
        ep <- e$win_stop[[k]]
        mc <- e$win_mc[[k]]
        spf <- sample_uniform_peptidoform_in_window(e$sequence, sp, ep, maps$aa_to_types)
        row <- data.frame(
          Accession = e$accession,
          Peptide = spf$Peptide,
          Start = spf$Start,
          Stop = spf$Stop,
          MC = mc,
          stringsAsFactors = FALSE
        )
        row$PTMPos <- spf$PTMPos
        row$PTMType <- spf$PTMType
        out[[out_ptr]] <- row
        out_ptr <- out_ptr + 1L
      }
    } else {
      # Fallback: uniform windows (e.g., no PTMs case or missing win_count)
      counts <- lengths(e$valid_mc)
      if (length(counts) == 0) next
      for (jj in seq_len(n_draws)) {
        s <- sample.int(length(counts), 1L, prob = counts)
        mc <- sample(e$valid_mc[[s]], 1L)
        sp <- e$starts[s]
        ep <- e$stops[s + mc]
        spf <- sample_uniform_peptidoform_in_window(e$sequence, sp, ep, maps$aa_to_types)
        row <- data.frame(
          Accession = e$accession,
          Peptide = spf$Peptide,
          Start = spf$Start,
          Stop = spf$Stop,
          MC = mc,
          stringsAsFactors = FALSE
        )
        row$PTMPos <- spf$PTMPos
        row$PTMType <- spf$PTMType
        out[[out_ptr]] <- row
        out_ptr <- out_ptr + 1L
      }
    }
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
  maps <- build_aa_maps(parameters)

  idx <- vector("list", nrow(df))
  dropped_no_windows <- 0L
  windows_per_protein <- integer(nrow(df))
  pep2prot_env <- new.env(parent = emptyenv())

  # Parallelization controls via environment (default serial)
  cores_env <- suppressWarnings(as.integer(Sys.getenv("PM_CORES", "1")))
  if (is.na(cores_env) || cores_env < 1) cores_env <- 1L
  cluster_env <- toupper(Sys.getenv("PM_CLUSTER", "PSOCK"))

  if (cores_env > 1L) {
    message(" + Precomputing per-protein windows and counts: parallel (", cluster_env, ", cores=", cores_env, ") over ", nrow(df), " proteins …")
    cl <- parallel::makeCluster(cores_env, type = cluster_env)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    # Export small, read-only objects
    aa_count <- maps$aa_to_count
    pep_min <- parameters$PepMinLength
    pep_max <- parameters$PepMaxLength
    max_mc  <- parameters$MaxNumMissedCleavages
    parallel::clusterExport(cl, varlist = c("aa_count", "cre", "pep_min", "pep_max", "max_mc"), envir = environment())
    # Build tasks
    tasks <- lapply(seq_len(nrow(df)), function(i) list(acc = df$Accession[i], seq = df$Sequence[i]))
    parts <- parallel::parLapply(cl, tasks, function(t) {
      # Per protein worker
      seqi <- t$seq
      sites <- stringi::stri_locate_all_regex(seqi, cre)[[1]]
      cs <- if (is.null(sites) || is.na(sites[1,1])) integer(0) else as.integer(sites[,1])
      if (length(cs) == 0) {
        starts <- 1L; stops <- nchar(seqi)
      } else {
        starts <- c(1L, cs + 1L)
        stops  <- c(cs, nchar(seqi))
      }
      S <- length(starts)
      if (S == 0) return(list(drop = TRUE))
      # Vectorized enumeration of windows using interval searches
      target_min <- starts + (pep_min - 1L)
      target_max <- starts + (pep_max - 1L)
      e_min <- pmax(seq_len(S), findInterval(target_min - 1L, stops) + 1L)
      e_max <- pmin(seq_len(S) + max_mc, findInterval(target_max, stops))
      counts <- e_max - e_min + 1L
      counts[counts < 0L] <- 0L
      n_valid <- sum(counts)
      if (n_valid == 0L) return(list(drop = TRUE))
      # Expand to window vectors
      s_rep <- rep.int(seq_len(S), counts)
      end_idx <- rep(e_min, counts) + sequence(counts) - 1L
      st_pos <- starts[s_rep]
      en_pos <- stops[end_idx]
      mc_v <- end_idx - s_rep
      # Reconstruct valid_mc list (length S) for compatibility
      valid_mc <- vector("list", S)
      if (length(mc_v) > 0) {
        split_list <- split(mc_v, s_rep)
        idx_names <- as.integer(names(split_list))
        valid_mc[idx_names] <- split_list
      }
      # ai_counts
      aa_chars <- strsplit(seqi, "", fixed = TRUE)[[1]]
      map_idx <- match(aa_chars, LETTERS)
      aic <- as.integer(aa_count[map_idx]); aic[is.na(aic)] <- 0L
      # per-window counts via log-prefix
      log_prefix <- c(0, cumsum(log1p(aic)))
      log_cnt <- log_prefix[en_pos + 1L] - log_prefix[st_pos]
      cnt_vec <- exp(log_cnt)
      list(
        drop = FALSE,
        accession = t$acc,
        sequence  = seqi,
        starts = starts, stops = stops, valid_mc = valid_mc, n_valid = n_valid,
        ai_counts = aic,
        win_start = st_pos, win_stop = en_pos, win_mc = as.integer(mc_v), win_count = as.numeric(cnt_vec),
        pf_total = sum(cnt_vec)
      )
    })
    # Assemble results
    for (i in seq_along(parts)) {
      p <- parts[[i]]
      if (is.null(p) || isTRUE(p$drop)) {
        idx[[i]] <- NULL
        dropped_no_windows <- dropped_no_windows + 1L
        windows_per_protein[i] <- 0L
        next
      }
      windows_per_protein[i] <- p$n_valid
      idx[[i]] <- list(
        accession = p$accession,
        sequence  = p$sequence,
        starts    = p$starts,
        stops     = p$stops,
        valid_mc  = p$valid_mc,
        n_valid   = p$n_valid,
        ai_counts = p$ai_counts,
        pf_total  = p$pf_total,
        win_start = p$win_start,
        win_stop  = p$win_stop,
        win_mc    = p$win_mc,
        win_count = p$win_count
      )
      # Build peptide -> protein map (I->L normalized)
      if (length(p$win_start) > 0) {
        pep <- substring(p$sequence, p$win_start, p$win_stop)
        build_pep2prot_map(normalize_IL(pep), p$accession, pep2prot_env)
      }
    }
  } else {
    # Serial path
    for (i in seq_len(nrow(df))) {
      seqi <- df$Sequence[i]
      sites <- cleavage_sites(seqi, cre)
      seg   <- segments_from_sites(seqi, sites)
      # Vectorized enumeration
      S <- length(seg$starts)
      if (S == 0) {
        idx[[i]] <- NULL
        dropped_no_windows <- dropped_no_windows + 1L
        next
      }
      target_min <- seg$starts + (parameters$PepMinLength - 1L)
      target_max <- seg$starts + (parameters$PepMaxLength - 1L)
      e_min <- pmax(seq_len(S), findInterval(target_min - 1L, seg$stops) + 1L)
      e_max <- pmin(seq_len(S) + parameters$MaxNumMissedCleavages, findInterval(target_max, seg$stops))
      counts <- e_max - e_min + 1L
      counts[counts < 0L] <- 0L
      n_valid <- sum(counts)
      windows_per_protein[i] <- n_valid
      if (n_valid == 0L) {
        idx[[i]] <- NULL
        dropped_no_windows <- dropped_no_windows + 1L
        next
      }

      # Precompute per-position ai counts (#PTM types per residue)
      aa_chars <- strsplit(seqi, "", fixed = TRUE)[[1]]
      ai_counts <- as.integer(maps$aa_to_count[aa_chars])
      # Sum of peptidoform counts across all valid peptide windows in this protein
      pf_total <- 0
      win_start <- integer(0)
      win_stop  <- integer(0)
      win_mc    <- integer(0)
      win_count <- numeric(0)
      # Expand to window vectors
      s_rep <- rep.int(seq_len(S), counts)
      end_idx <- rep(e_min, counts) + sequence(counts) - 1L
      st_pos <- seg$starts[s_rep]
      en_pos <- seg$stops[end_idx]
      MC <- end_idx - s_rep
      # Use log-prefix to avoid overflow, then exponentiate per window
      ai_counts[is.na(ai_counts)] <- 0L
      log_prefix <- c(0, cumsum(log1p(ai_counts)))
      log_cnt <- log_prefix[en_pos + 1L] - log_prefix[st_pos]
      cnt_vec <- exp(log_cnt)
      pf_total <- sum(cnt_vec)
      win_start <- st_pos
      win_stop  <- en_pos
      win_mc    <- MC
      win_count <- as.numeric(cnt_vec)

      idx[[i]] <- list(
        accession = df$Accession[i],
        sequence  = seqi,
        starts    = seg$starts,
        stops     = seg$stops,
        valid_mc  = split(win_mc, s_rep),
        n_valid   = n_valid,
        ai_counts = ai_counts,
        pf_total  = pf_total,
        win_start = win_start,
        win_stop  = win_stop,
        win_mc    = win_mc,
        win_count = win_count
      )

      # Build peptide -> protein map (I->L normalized)
      if (length(win_start) > 0) {
        pep <- substring(seqi, win_start, win_stop)
        build_pep2prot_map(normalize_IL(pep), df$Accession[i], pep2prot_env)
      }
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
  # Optional summary of peptidoform totals per protein (sum over windows of prod(1+a_i))
  if (length(idx) > 0 && !is.null(idx[[1]]$pf_total)) {
    pf_totals <- vapply(idx, function(x) as.numeric(x$pf_total), numeric(1))
    cat(" + Peptidoform totals per protein (non-zero windows): min=", round(min(pf_totals),2),
        " median=", round(stats::median(pf_totals),2),
        " mean=", round(mean(pf_totals),2),
        " max=", round(max(pf_totals),2), "\n", sep = "")
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

# Report total index build time (explicit)
if (!is.null(res$build_time_sec)) {
  cat(sprintf("\n#INDEX BUILD - Done (%.3f sec)\n", as.numeric(res$build_time_sec)))
}

if (!is.na(cli$outRDS)) {
  saveRDS(res, file = cli$outRDS)
  cat("\nSaved index to:", cli$outRDS, "\n")
}

# Optional: sample a random set of peptidoforms uniformly across the full space
if (!is.na(cli$subsetN) && cli$subsetN > 0) {
  cat("\n#PEPTIDOFORMS SUBSET - Start\n")
  # Diagnostics to help debug issues before sampling
  n_prot <- length(res$proteins)
  n_valid_vec <- vapply(res$proteins, function(e) if (is.null(e)) 0L else e$n_valid, integer(1))
  pf_vec <- vapply(res$proteins, function(e) if (is.null(e) || is.null(e$pf_total)) 0 else as.numeric(e$pf_total), numeric(1))
  invalid_win <- vapply(res$proteins, function(e) {
    wc <- if (is.null(e)) NULL else e$win_count
    if (is.null(wc) || length(wc) == 0) return(FALSE)
    any(!is.finite(wc)) || any(is.na(wc))
  }, logical(1))
  ai_sum <- sum(vapply(res$proteins, function(e) if (is.null(e) || is.null(e$ai_counts)) 0L else sum(e$ai_counts, na.rm = TRUE), numeric(1)))
  cat(sprintf(" + Subset diagnostics: proteins=%d, total_windows=%d, zero_window_proteins=%d\n",
              n_prot, sum(n_valid_vec), sum(n_valid_vec == 0)))
  cat(sprintf(" + pf_total: sum=%.0f, zero_pf_total_proteins=%d, invalid_window_counts=%d\n",
              sum(pf_vec), sum(pf_vec == 0), sum(invalid_win)))
  cat(sprintf(" + PTM configuration: total sum of a_i across proteome=%d (%s)\n",
              as.integer(ai_sum), if (ai_sum == 0) "no PTMs (uniform over windows)" else "PTMs present"))
  t_sub0 <- proc.time()[[3]]
  subset <- sample_peptidoforms(res, params, cli$subsetN)
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
  t_sub1 <- proc.time()[[3]]
  cat(sprintf("#PEPTIDOFORMS SUBSET - Finish (%.3f sec)\n", t_sub1 - t_sub0))
}
