################################################################################
#                    PEPTIDE DIGESTION OF PROTEOFORM TABLE                     #
################################################################################


# Internal helpers for enzyme cleavage and indexing (package-internal)
# These utilities can be reused by digestion and potential indexing/caching.
# They are intentionally lightweight and keep behavior consistent with fastDigest.

#' @keywords internal
enzymeRegex <- function(enzyme) {
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

#' @keywords internal
cleavageSites <- function(sequence, cre) {
  loc <- stringi::stri_locate_all_regex(sequence, cre)[[1]]
  if (is.null(loc) || is.na(loc[1, 1])) integer(0) else as.integer(loc[, 1])
}

#' @keywords internal
segmentsFromSites <- function(sequence, sites) {
  if (length(sites) == 0) {
    list(starts = 1L, stops = nchar(sequence))
  } else {
    list(starts = c(1L, sites + 1L), stops = c(sites, nchar(sequence)))
  }
}

#' @keywords internal
enumerateValidWindows <- function(starts, stops, pep_min, pep_max, max_mc) {
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

#' Build a lightweight search index from in-memory sequences
#'
#' @param proteins data.frame with columns Accession and Sequence
#' @param parameters list with Enzyme, PepMinLength, PepMaxLength, MaxNumMissedCleavages
#' @param includePep2Prot logical, whether to compute shared peptide map
#' @return list with proteins, weights, optional pep2prot, params, build_time_sec
#' @keywords internal
buildSearchIndexFromSequences <- function(proteins, parameters) {
  t0 <- proc.time()[[3]]
  cre <- enzymeRegex(parameters$Enzyme)

  idx <- vector("list", nrow(proteins))
  dropped_no_windows <- 0L
  windows_per_protein <- integer(nrow(proteins))
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

  # Optional parallelization using parameters$Cores and parameters$ClusterType
  cores <- tryCatch(parameters$Cores, error = function(e) NULL)
  cl_type <- toupper(tryCatch(parameters$ClusterType, error = function(e) ""))

  # Build AA count map once for worker export
  # (mirror of buildAAMapsDigest but only counts are needed here)
  local_aa_count <- {
    default_types <- c("ph", "ox", "ac", "me", "de")
    default_residues <- list(
      ph = c("S", "T", "Y"),
      ox = c("M", "W", "C", "Y"),
      ac = c("K"),
      me = c("K", "R"),
      de = c("D", "E")
    )
    ptm_types <- parameters$PTMTypes
    if (is.null(ptm_types) || length(ptm_types) == 0 || all(is.na(ptm_types))) ptm_types <- list(default_types)
    if (is.list(ptm_types) && length(ptm_types) == 1) ptm_types <- ptm_types[[1]]
    ptm_types <- ptm_types[!is.na(ptm_types)]
    modres <- parameters$ModifiableResidues
    if (is.null(modres) || length(modres) == 0 || all(is.na(modres))) {
      modres_map <- default_residues
    } else if (is.list(modres) && length(modres) == 1 && is.list(modres[[1]])) {
      modres_map <- modres[[1]]
    } else if (is.list(modres)) {
      modres_map <- modres
    } else {
      modres_map <- default_residues
    }
    aa_to_types <- setNames(vector("list", length = 26L), LETTERS)
    for (ptm in ptm_types) {
      aa <- modres_map[[ptm]]
      if (!is.null(aa) && length(aa) > 0) for (a in aa) aa_to_types[[a]] <- unique(c(aa_to_types[[a]], ptm))
    }
    aa_to_count <- setNames(integer(26L), LETTERS)
    for (a in names(aa_to_types)) aa_to_count[[a]] <- length(aa_to_types[[a]])
    aa_to_count
  }
  aa_count <- local_aa_count

  if (!is.null(cores) && is.numeric(cores) && cores > 1L && cl_type %in% c("PSOCK","FORK")) {
    cl <- parallel::makeCluster(min(cores, parallel::detectCores()), type = cl_type)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    pep_min <- parameters$PepMinLength
    pep_max <- parameters$PepMaxLength
    max_mc  <- parameters$MaxNumMissedCleavages
    parallel::clusterExport(cl, varlist = c("cre","pep_min","pep_max","max_mc","aa_count"), envir = environment())
    tasks <- lapply(seq_len(nrow(proteins)), function(i) list(acc = proteins$Accession[i], seq = proteins$Sequence[i]))
    parts <- parallel::parLapply(cl, tasks, function(t) {
      # Cleavage → segments
      loc <- stringi::stri_locate_all_regex(t$seq, cre)[[1]]
      cs <- if (is.null(loc) || is.na(loc[1,1])) integer(0) else as.integer(loc[,1])
      if (length(cs) == 0) {
        starts <- 1L; stops <- nchar(t$seq)
      } else {
        starts <- c(1L, cs + 1L); stops <- c(cs, nchar(t$seq))
      }
      S <- length(starts); if (S == 0) return(NULL)
      # Vectorized windows
      target_min <- starts + (pep_min - 1L)
      target_max <- starts + (pep_max - 1L)
      e_min <- pmax(seq_len(S), findInterval(target_min - 1L, stops) + 1L)
      e_max <- pmin(seq_len(S) + max_mc, findInterval(target_max, stops))
      counts <- e_max - e_min + 1L; counts[counts < 0L] <- 0L
      n_valid <- sum(counts); if (n_valid == 0L) return(NULL)
      s_rep <- rep.int(seq_len(S), counts)
      end_idx <- rep(e_min, counts) + sequence(counts) - 1L
      st_pos <- starts[s_rep]; en_pos <- stops[end_idx]; mc_v <- end_idx - s_rep
      # ai_counts and per-window counts
      aic <- as.integer(aa_count[match(strsplit(t$seq, "", fixed = TRUE)[[1]], LETTERS)])
      aic[is.na(aic)] <- 0L
      lp <- c(0, cumsum(log1p(aic)))
      cnt_vec <- exp(lp[en_pos + 1L] - lp[st_pos])
      # Reconstruct valid_mc list
      valid_mc <- vector("list", S)
      if (length(mc_v) > 0) {
        split_list <- split(mc_v, s_rep)
        idx_names <- as.integer(names(split_list))
        valid_mc[idx_names] <- split_list
      }
      list(
        accession = t$acc,
        sequence  = t$seq,
        starts = starts, stops = stops, n_valid = n_valid,
        ai_counts = aic,
        win_start = st_pos, win_stop = en_pos, win_mc = as.integer(mc_v), win_count = as.numeric(cnt_vec),
        pf_total = sum(cnt_vec)
      )
    })
    for (i in seq_along(parts)) {
      p <- parts[[i]]
      if (is.null(p)) {
        idx[[i]] <- NULL; dropped_no_windows <- dropped_no_windows + 1L; windows_per_protein[i] <- 0L; next
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
      if (length(p$win_start) > 0) {
        pep <- substring(p$sequence, p$win_start, p$win_stop)
        add_map(gsub("I", "L", pep, perl = TRUE), p$accession)
      }
    }
  } else for (i in seq_len(nrow(proteins))) {
    seqi <- proteins$Sequence[i]
    sites <- cleavageSites(seqi, cre)
    seg   <- segmentsFromSites(seqi, sites)

    # Vectorized enumeration of valid windows
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

    # Expand to per-window vectors
    s_rep   <- rep.int(seq_len(S), counts)
    end_idx <- rep(e_min, counts) + sequence(counts) - 1L
    win_start <- seg$starts[s_rep]
    win_stop  <- seg$stops[end_idx]
    win_mc    <- end_idx - s_rep

    # Per-residue PTM counts and per-window peptidoform counts via log-prefix
    aa_chars <- strsplit(seqi, "", fixed = TRUE)[[1]]
    # Build per-AA PTM count map (use existing maps below)
    # Reuse aa_maps built later; compute locally here
    # Build small map on the fly: default handled by log1p(0)
    # For speed, derive counts via string match with ModifiableResidues would be heavier; use existing mapping below
    # We rebuild a minimal map here:
    # Simpler: get from parameters via buildAAMapsDigest in this scope
    # However, aa_maps is created later; reusing buildAAMapsDigest now
    aac_map <- tryCatch(buildAAMapsDigest(parameters)$aa_to_count, error = function(e) NULL)
    if (is.null(aac_map)) {
      aac_map <- setNames(integer(26L), LETTERS)
    }
    ai_counts <- as.integer(aac_map[match(aa_chars, LETTERS)])
    ai_counts[is.na(ai_counts)] <- 0L
    log_prefix <- c(0, cumsum(log1p(ai_counts)))
    log_cnt <- log_prefix[win_stop + 1L] - log_prefix[win_start]
    win_count <- as.numeric(exp(log_cnt))
    pf_total <- sum(win_count)

    # valid_mc reconstructed for compatibility
    valid_mc <- vector("list", S)
    if (length(win_mc) > 0) {
      split_list <- split(win_mc, s_rep)
      idx_names <- as.integer(names(split_list))
      valid_mc[idx_names] <- split_list
    }

    idx[[i]] <- list(
      accession = proteins$Accession[i],
      sequence  = seqi,
      starts    = seg$starts,
      stops     = seg$stops,
      valid_mc  = valid_mc,
      n_valid   = n_valid,
      ai_counts = ai_counts,
      pf_total  = pf_total,
      win_start = win_start,
      win_stop  = win_stop,
      win_mc    = win_mc,
      win_count = win_count
    )

    # Peptide -> protein map (I->L normalized) from window vectors
    if (length(win_start) > 0) {
      pep <- substring(seqi, win_start, win_stop)
      pep_norm <- gsub("I", "L", pep, perl = TRUE)
      add_map(pep_norm, proteins$Accession[i])
    }
  }

  keep_mask <- !vapply(idx, is.null, logical(1))
  idx <- idx[keep_mask]
  # Accession -> index lookup for O(1) access during digestion
  acc2idx <- integer(0)
  if (length(idx) > 0) {
    accs <- vapply(idx, function(x) x$accession, character(1))
    acc2idx <- stats::setNames(seq_along(accs), accs)
  }
  raw_w <- if (length(idx) > 0) vapply(idx, function(x) x$n_valid, integer(1)) else integer(0)
  weights <- if (length(raw_w) == 0) {
    numeric(0)
  } else if (sum(raw_w) > 0) {
    raw_w / sum(raw_w)
  } else {
    rep(1 / length(raw_w), length(raw_w))
  }

  # Precompute peptidoform window probabilities for fast donor sampling in MS stage
  buildAAMapsDigest <- function(parameters) {
    default_types <- c("ph", "ox", "ac", "me", "de")
    default_residues <- list(
      ph = c("S", "T", "Y"),
      ox = c("M", "W", "C", "Y"),
      ac = c("K"),
      me = c("K", "R"),
      de = c("D", "E")
    )
    ptm_types <- parameters$PTMTypes
    if (is.null(ptm_types) || length(ptm_types) == 0 || all(is.na(ptm_types))) ptm_types <- list(default_types)
    if (is.list(ptm_types) && length(ptm_types) == 1) ptm_types <- ptm_types[[1]]
    ptm_types <- ptm_types[!is.na(ptm_types)]
    modres <- parameters$ModifiableResidues
    if (is.null(modres) || length(modres) == 0 || all(is.na(modres))) {
      modres_map <- default_residues
    } else if (is.list(modres) && length(modres) == 1 && is.list(modres[[1]])) {
      modres_map <- modres[[1]]
    } else if (is.list(modres)) {
      modres_map <- modres
    } else {
      modres_map <- default_residues
    }
    aa_to_types <- setNames(vector("list", length = 26L), LETTERS)
    for (ptm in ptm_types) {
      aa <- modres_map[[ptm]]
      if (!is.null(aa) && length(aa) > 0) for (a in aa) aa_to_types[[a]] <- unique(c(aa_to_types[[a]], ptm))
    }
    aa_to_count <- setNames(integer(26L), LETTERS)
    for (a in names(aa_to_types)) aa_to_count[[a]] <- length(aa_to_types[[a]])
    list(aa_to_types = aa_to_types, aa_to_count = aa_to_count)
  }
  aa_maps <- buildAAMapsDigest(parameters)
  t1 <- proc.time()[[3]]

  pep2prot <- list()
  if (length(ls(envir = pep2prot_env, all.names = TRUE)) > 0) {
    keys <- ls(envir = pep2prot_env, all.names = TRUE)
    vals <- lapply(keys, function(k) unique(get(k, envir = pep2prot_env, inherits = FALSE)))
    keep <- vapply(vals, function(v) length(v) > 1, logical(1))
    if (any(keep)) {
      pep2prot <- stats::setNames(lapply(vals[keep], function(x) sort(x)), keys[keep])
    }
  }

  list(
    proteins = idx,
    acc2idx = acc2idx,
    aa_to_types = aa_maps$aa_to_types,
    weights = weights,
    pep2prot = pep2prot,
    params = parameters,
    build_time_sec = t1 - t0,
    dropped_no_windows = dropped_no_windows,
    windows_per_protein = windows_per_protein
  )
}

#' Build a search index directly from FASTA on disk (full proteome)
#'
#' Uses the same logic as BuildSearchIndex CLI but stays in-package.
#' @param parameters list with PathToFasta, Enzyme, PepMinLength, PepMaxLength, MaxNumMissedCleavages
#' @return index list as in buildSearchIndexFromSequences
#' @keywords internal
buildSearchIndexFromFasta <- function(parameters) {
  t_start <- proc.time()[[3]]
  message("\n#PROTEOME INDEX - Start\n")
  message(" + FASTA: ", parameters$PathToFasta)
  message(" + Enzyme: ", parameters$Enzyme,
          ", len ", parameters$PepMinLength, "-", parameters$PepMaxLength,
          ", maxMC=", parameters$MaxNumMissedCleavages)

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
    keep <- !Reduce(`|`, lapply(badAA, function(x) grepl(x, df$Sequence, fixed = TRUE)))
    df <- df[keep, , drop = FALSE]
  }
  df <- df[!duplicated(df$Accession), , drop = FALSE]
  message(" + Proteins after filtering: ", nrow(df))

  idx <- buildSearchIndexFromSequences(df, parameters)

  total_windows <- if (length(idx$proteins) > 0) sum(vapply(idx$proteins, function(x) x$n_valid, integer(1))) else 0L
  wp <- idx$windows_per_protein[idx$windows_per_protein > 0]
  med_wp <- if (length(wp) > 0) stats::median(wp) else NA
  mean_wp <- if (length(wp) > 0) round(mean(wp), 2) else NA
  min_wp <- if (length(wp) > 0) min(wp) else NA
  max_wp <- if (length(wp) > 0) max(wp) else NA
  n_shared <- length(idx$pep2prot)
  t_total <- proc.time()[[3]] - t_start

  message(" + Indexed proteins: ", length(idx$proteins))
  message(" + Dropped (no windows): ", idx$dropped_no_windows)
  message(" + Total valid windows: ", total_windows)
  message(" + Windows/protein (nz) min/med/mean/max: ", min_wp, "/", med_wp, "/", mean_wp, "/", max_wp)
  message(" + Shared peptides (>=2 proteins): ", n_shared)
  message(sprintf(" + Time: %.3f sec", t_total))
  message("#PROTEOME INDEX - Finish\n")

  idx
}



#' Perform enzymatic digestion on a single protein sequence
#'
#' This function simulates enzymatic digestion on a single protein sequence, supporting multiple
#' cleavage rules for enzymes such as trypsin, chymotrypsin, pepsin, lysC, and argC. The function
#' returns peptides filtered by length and mass, including calculations for charges +1, +2, and +3.
#'
#' @param sequence A character string representing the protein sequence to be digested.
#' @param enzyme A character string specifying the enzyme to use for digestion. Default is "trypsin".
#' @param missed An integer specifying the maximum number of missed cleavages. Default is 0.
#' @param length.max An integer specifying the maximum peptide length. Default is NA.
#' @param length.min An integer specifying the minimum peptide length. Default is NA.
#' @param mass.max A numeric value specifying the maximum peptide mass. Default is NA.
#' @param mass.min A numeric value specifying the minimum peptide mass. Default is NA.
#'
#' @return A data frame with the digested peptides, including their mass/charge (M/Z) ratios for charges +1, +2, and +3.
#' If no peptides meet the criteria, the function returns \code{NULL}.
#'
#' @importFrom dplyr case_when bind_rows
#' @importFrom stringi stri_locate_all_regex
#' @importFrom crayon red
#' @keywords internal
fastDigest <- function(proteoform, parameters, searchIndex) {
  acc <- as.character(proteoform$Accession)
  if (length(acc) != 1L || is.na(acc) || !nzchar(acc)) {
    message("[fastDigest] Skipping proteoform with invalid accession: ", acc)
    return(NULL)
  }
  # Find index entry for this accession (index is assumed to be present)
  e <- NULL
  if (!is.null(searchIndex$acc2idx)) {
    idx <- searchIndex$acc2idx[[acc]]
    if (!is.null(idx)) {
      e <- searchIndex$proteins[[ idx ]]
    } else {
      message("[fastDigest] Accession not found in index: ", acc)
      return(NULL)
    }
  } else {
    # Fallback linear scan (should be rare)
    idx_match <- which(vapply(searchIndex$proteins, function(x) identical(x$accession, acc), logical(1)))
    if (length(idx_match) == 0) return(NULL)
    e <- searchIndex$proteins[[idx_match[1]]]
  }

  # Use precomputed window vectors from the index
  if (is.null(e$win_start) || length(e$win_start) == 0) return(NULL)
  peptides <- data.frame(
    Peptide = substring(e$sequence, e$win_start, e$win_stop),
    Start = e$win_start,
    Stop  = e$win_stop,
    MC    = e$win_mc,
    stringsAsFactors = FALSE
  )

  if (nrow(peptides) == 0) return(NULL)

  AAs <- strsplit(peptides$Peptide, split = "")
  AA.mass <- c(
    "A" = 71.03711, "R" = 156.10111, "N" = 114.04293, "D" = 115.02694, "C" = 103.00919,
    "E" = 129.04259, "Q" = 128.05858, "G" = 57.02146, "H" = 137.05891, "I" = 113.08406,
    "L" = 113.08406, "K" = 128.09496, "M" = 131.04049, "F" = 147.06841, "P" = 97.05276,
    "S" = 87.03203, "T" = 101.04768, "W" = 186.07931, "Y" = 163.06333, "V" = 99.06841
  )
  peptide.mass <- sapply(AAs, function(x) sum(AA.mass[x], 18.01528))
  peptides$MZ1 <- peptide.mass + 1.007276466
  peptides$MZ2 <- (peptide.mass + (1.007276466 * 2)) / 2
  peptides$MZ3 <- (peptide.mass + (1.007276466 * 3)) / 3

  rownames(peptides) <- if (nrow(peptides) > 0) 1:nrow(peptides) else NULL
  peptides
}
#####################

#' Perform enzymatic digestion on a set of proteoforms
#'
#' This function performs enzymatic digestion on a set of proteoforms, mapping modification sites
#' on peptide sequences and adding mass shifts per peptide based on modifications. The peptide
#' abundance is set based on the parental proteoform.
#'
#' @param proteoform A data frame containing proteoform sequences and associated data.
#' @param parameters A list containing various parameters for the digestion process,
#' including enzyme type, maximum number of missed cleavages, peptide length limits,
#' and modification masses.
#'
#' @return A data frame containing the digested peptides with mapped modifications,
#' mass shifts, and associated abundances. If no peptides are generated, the function
#' returns \code{NULL}.
#'
#' @importFrom dplyr bind_cols
#' @importFrom stats aggregate
#' @keywords internal
proteoformDigestion <- function(proteoform, parameters, searchIndex = NULL) {
  peptides <- fastDigest(proteoform = proteoform, parameters = parameters, searchIndex = searchIndex)

  if (!is.null(peptides)) {
    peptides$Accession <- proteoform$Accession
    peptides$Proteoform_ID <- proteoform$Proteoform_ID

    peptides$PTMPos <- vector(mode = "list", length = nrow(peptides))
    peptides$PTMType <- vector(mode = "list", length = nrow(peptides))
    peptides$Regulation_Amplitude <- proteoform$Regulation_Amplitude
    peptides$Regulation_Pattern <- proteoform$Regulation_Pattern

    # Map modification sites on peptides.
    if (!is.null(proteoform$PTMPos[[1]])) {
      proteoform.position <- unlist(proteoform$PTMPos)
      proteoform.type <- unlist(proteoform$PTMType)
      pep.indices <- lapply(proteoform.position, function(x) which(x >= peptides$Start & x <= peptides$Stop))

      if (sum(lengths(pep.indices)) != 0) {
        pep.position <- lapply(1:length(pep.indices), function(x) sapply(pep.indices[[x]], function(y) proteoform.position[x] - peptides$Start[y] + 1))
        pep.type <- lapply(1:length(pep.indices), function(x) rep(proteoform.type[x], length(pep.indices[[x]])))

        to.aggregate <- data.frame(unlist(pep.indices), unlist(pep.position), unlist(pep.type), stringsAsFactors = F)
        to.aggregate <- stats::aggregate(to.aggregate[, 2:3], by = list(to.aggregate[, 1]), FUN = list)

        # Calculate and add the mass addition due to modifications per modified peptide.
        modification.mass <- parameters$PTMTypesMass[[1]]
        names(modification.mass) <- parameters$PTMTypes[[1]]
        to.aggregate$mass_shift <- sapply(to.aggregate[, 3], function(x) sum(unlist(modification.mass[x]), na.rm = T))

        peptides[to.aggregate[, 1], c("PTMPos", "PTMType")] <- to.aggregate[, 2:3]
        peptides[to.aggregate[, 1], c("MZ1", "MZ2", "MZ3")] <- peptides[to.aggregate[, 1], c("MZ1", "MZ2", "MZ3")] + as.numeric(to.aggregate[, 4]) %*% t(c(1, 0.5, 1 / 3))
      }
    }

    # Add proteoform abundance to all peptides.
    peptides.abundance <- as.data.frame(matrix(NA, ncol = length(parameters$QuantColnames), nrow = nrow(peptides)))
    colnames(peptides.abundance) <- parameters$QuantColnames
    peptides.abundance[1:nrow(peptides.abundance), parameters$QuantColnames] <- proteoform[parameters$QuantColnames]

    # Bind everything.
    peptides <- dplyr::bind_cols(peptides, peptides.abundance)
  }

  return(peptides)
}
#####################

#' Perform proteoform digestion on a set of proteoforms with optional parallel computing
#'
#' This function wraps the \code{proteoformDigestion} function to perform enzymatic digestion
#' on a set of proteoforms. It supports parallel computing and allows sampling of peptides
#' based on a distribution derived from \code{PropMissedCleavages} according to their missed
#' cleavages (MC).
#'
#' @param proteoforms A data frame containing the proteoform sequences and associated data.
#' @param parameters A list of parameters including enzyme type, number of cores for parallel
#' computing, maximum number of missed cleavages, peptide length limits, and proportion of missed cleavages.
#'
#' @return A data frame containing the digested peptides, with details on modifications,
#' mass shifts, and missed cleavages.
#'
#' @importFrom dplyr bind_rows
#' @importFrom parallel makeCluster detectCores setDefaultCluster clusterExport stopCluster parLapply
#' @importFrom scales rescale
#' @keywords internal
digestGroundTruth <- function(proteoforms, parameters, searchIndex = NULL) {
  message("\n#PROTEOFORM DIGESTION - Start\n")
  message(" + Digestion input:")
  message("  - A total number of ", nrow(proteoforms), " proteoforms, is proceed for proteolytic digestion.")
  message("  - Unmodified fraction contains ", sum(lengths(proteoforms$PTMType) == 0), " proteoforms and modified fraction ", sum(lengths(proteoforms$PTMType) != 0), " proteoforms.")
  message(
    "  - Cleavage will be performed by ", parameters$Enzyme, " with a maximum of ", parameters$MaxNumMissedCleavages,
    " miss-cleavages, to create peptides of length ", parameters$PepMinLength, " to ", parameters$PepMaxLength, " amino acids."
  )

  # Filter proteoforms whose parent protein has no valid peptide windows in the index (if provided)
  if (!is.null(searchIndex)) {
    acc2idx <- searchIndex$acc2idx
    accs <- as.character(proteoforms$Accession)
    accs[!nzchar(accs)] <- NA_character_
    idx_vec <- unname(acc2idx[accs])
    missing <- sum(is.na(idx_vec))
    valid_mask <- rep(FALSE, length(idx_vec))
    good <- which(!is.na(idx_vec))
    if (length(good)) {
      valid_mask[good] <- vapply(idx_vec[good], function(j) {
        e <- searchIndex$proteins[[j]]
        !is.null(e$win_start) && length(e$win_start) > 0
      }, logical(1))
    }
    num_discard <- sum(!valid_mask)
    if (num_discard > 0) {
      message("  - Discarding ", num_discard, " proteoforms without valid peptide windows in index (", missing, " missing accessions).")
      proteoforms <- proteoforms[valid_mask, , drop = FALSE]
    }
    if (nrow(proteoforms) == 0) {
      message("  - No proteoforms left to digest after filtering.\n#PROTEOFORM DIGESTION - Finish\n")
      return(NULL)
    }
  }

  # Prefer serial digestion unless FORK clusters are used (to avoid copying large index)
  if (!is.null(parameters$Cores) && parameters$Cores > 1 && identical(tolower(parameters$ClusterType), "fork")) {
    cores <- parameters$Cores

    if (parallel::detectCores() <= parameters$Cores) {
      cores <- parallel::detectCores() - 1
    }

    cluster <- parallel::makeCluster(cores, type = parameters$ClusterType)
    #on.exit(parallel::stopCluster(cluster))
    parallel::setDefaultCluster(cluster)
    parallel::clusterExport(cluster, c("proteoforms", "parameters", "proteoformDigestion", "fastDigest", "searchIndex"), envir = environment())
    peptides <- parallel::parLapply(cluster, 1:nrow(proteoforms), function(x) proteoformDigestion(proteoform = proteoforms[x, ], parameters = parameters, searchIndex = searchIndex))
    parallel::stopCluster(cluster)
  } else {
    peptides <- lapply(1:nrow(proteoforms), function(x) proteoformDigestion(proteoform = proteoforms[x, ], parameters = parameters, searchIndex = searchIndex))
  }

  message("  - All proteoforms are digested successfully!")
  # Discard proteoforms without valid peptide windows or missing from index
  num_discarded <- sum(vapply(peptides, is.null, logical(1)))
  if (num_discarded > 0) {
    message("  - Discarded ", num_discarded, " proteoforms without valid peptide windows (or not found in index).")
  }
  peptides <- peptides[!vapply(peptides, is.null, logical(1))]
  peptides <- dplyr::bind_rows(peptides)

  # Sample peptides per MC by size determined by PropMissedCleavages.
  if (parameters$MaxNumMissedCleavages > 0) {
    if (parameters$PropMissedCleavages > 0 & parameters$PropMissedCleavages < 1) {
      # set max number of missed cleavages for probability calculation to min 5
      max_misscleav <- ifelse(parameters$MaxNumMissedCleavages < 5, 5, parameters$MaxNumMissedCleavages)
      MC.proportions <- sapply(0:parameters$MaxNumMissedCleavages, function(x) choose(max_misscleav, x) *
                                 (1 - parameters$PropMissedCleavages)^(max_misscleav-x) *
                                 parameters$PropMissedCleavages^(x))
      # MC.proportions <- scales::rescale(x = MC.proportions, to = c(0, 1), from = c(0, max(MC.proportions, na.rm = T)))
      MC.proportions <- MC.proportions / max(MC.proportions)
      peptide.indices <- lapply(0:parameters$MaxNumMissedCleavages, function(x) which(peptides$MC == x))
      peptide.indices <- unlist(lapply(0:parameters$MaxNumMissedCleavages, function(x) {
        n_available <- length(peptide.indices[[x+1]])
        n_wanted <- floor(sum(peptides$MC == 0) * MC.proportions[x + 1])
        n_sample <- min(n_available, n_wanted)  # Clip
        if (n_sample > 0) sample(peptide.indices[[x+1]], size = n_sample, replace = FALSE) else NULL
      }))
      #      peptide.indices <- unlist(lapply(1:parameters$MaxNumMissedCleavages, function(x) sample(peptide.indices[[x]], size = floor(sum(peptides$MC == 0) * MC.proportions[x + 1]), replace = FALSE)))
      # peptide.indices <- sort(c(which(peptides$MC == 0), peptide.indices))
      peptides <- peptides[peptide.indices, ]
    } else if (parameters$PropMissedCleavages == 0) {
      peptides <- peptides[peptides$MC == 0, ]
    }
  }

  message(" + Digestion output:")
  message("  - A total number of ", nrow(peptides), " peptides is generated.")
  message("  - Unmodified fraction contains ", sum(lengths(peptides$PTMType) == 0), " peptides and modified fraction ", sum(lengths(peptides$PTMType) != 0), " peptides.")
  message(
    "  - The amount of peptides with ", paste0(0:parameters$MaxNumMissedCleavages, collapse = ", "), " miss-cleavages is ",
    paste0(sapply(0:parameters$MaxNumMissedCleavages, function(x) sum(peptides$MC == x)), collapse = ", "), " respectively.\n"
  )

  message("#PROTEOFORM DIGESTION - Finish\n")

  return(peptides)
}
#####################


#' Calculate the detectability of a peptide sequence
#'
#' This function calculates the detectability of a peptide sequence based on the amino acid
#' using the PeptideRanger package. It returns a numeric vector representing the detectability
#' of the peptides. The prediction of the detectability is based on the amino acid composition
#' and does not take into account post-translational modifications.
#'
#' @param peptides A character vector containing the peptide sequences.
#' @param parameters A list of parameters including the number of cores for parallel computing.
#'
#' @return A numeric vector representing the detectability of the peptides.
#'
#' @importFrom PeptideRanger peptide_predictions
#' @keywords internal
addDetectability <- function(peptides, parameters) {

  # get unique peptide list and be able to map back
  unique_peptides <- unique(peptides)
  peptide_map <- match(peptides, unique_peptides)

  RFScores <- NULL
  if (!is.null(parameters$Cores)) {
    cores <- parameters$Cores

    if (parallel::detectCores() <= parameters$Cores) {
      cores <- parallel::detectCores() - 1
    }

    cluster <- parallel::makeCluster(cores, type = parameters$ClusterType)
    parallel::setDefaultCluster(cluster)

    # Ensure that the necessary package is loaded on each worker
    parallel::clusterEvalQ(cluster, library(PeptideRanger))

    # Split the data into chunks of 100 peptides
    peptide_chunks <- split(unlist(unique_peptides) , ceiling(seq_along(unlist(unique_peptides))/100))

    # Run the predictions in parallel
    RFScores <- parallel::parLapply(cluster, peptide_chunks, function(subset) {
      PeptideRanger::peptide_predictions(unlist(subset), PeptideRanger::RFmodel_ProteomicsDB)
    })

    # Combine the results into a single list or data frame
    RFScores <- do.call(rbind, RFScores)
    parallel::stopCluster(cluster)
  } else {
    RFScores<- PeptideRanger::peptide_predictions(unlist(unique_peptides), PeptideRanger::RFmodel_ProteomicsDB)
  }

  # Map back to the original peptides
  RFScores <- RFScores[peptide_map,]$RF_score
  return(RFScores)
}
#####################


#' Summarize digested peptide products
#'
#' This function groups peptides by unique identifiers, summarizes their abundance,
#' and optionally removes a percentage of the least abundant peptides. It creates
#' unique peptide IDs, substitutes isoleucine with leucine, and aggregates peptides
#' based on various characteristics.
#'
#' @param peptides A data frame containing the digested peptides and associated data.
#' @param parameters A list of parameters including QuantColnames and LeastAbundantLoss.
#'
#' @return A data frame containing the summarized peptides, with the following structure:
#' \describe{
#'   \item{Sequence}{A character vector containing the unique peptide sequence after isoleucine substitution to leucine.}
#'   \item{Peptide}{A list of character vectors containing the peptides that are grouped based on Sequence, prior to isoleucine substitution.}
#'   \item{Start}{A list of integer vectors containing the starting positions of the peptides in the Peptide vectors on the protein sequence.}
#'   \item{Stop}{A list of integer vectors containing the ending positions of the peptides in the Peptide vectors on the protein sequence.}
#'   \item{MC}{A list of integer vectors containing the number of missed cleavages (MC) for the peptides in the Peptide vectors.}
#'   \item{MZ1}{A numeric vector representing the peptide mass for charge +1.}
#'   \item{MZ2}{A numeric vector representing the peptide mass for charge +2.}
#'   \item{MZ3}{A numeric vector representing the peptide mass for charge +3.}
#'   \item{Accession}{A list of character vectors containing the parental protein Accession of the peptides in the Peptide vectors.}
#'   \item{Proteoform_ID}{A list of integer vectors containing the unique proteoform identifiers of the Accession vectors.}
#'   \item{PTMPos}{A list of integer vectors containing the positions of the modifications on the peptides in the Peptide vectors.}
#'   \item{PTMType}{A list of character vectors containing the modification types of the modifications in the PTMPos vectors.}
#'   \item{Regulation_Amplitude}{A list of numeric vectors containing the regulation amplitudes of the proteoforms in the Accession vectors.}
#'   \item{Regulation_Pattern}{A list of numeric vectors containing the regulation patterns of the proteoforms in the Accession vectors.}
#'   \item{Quantitative Columns}{Numeric columns containing the abundances of the peptide group for each QuantColname. These columns are dynamically named based on the provided QuantColnames parameter.}
#' }
#'
#' @importFrom dplyr group_by summarise summarise_at vars inner_join select %>%
#' @keywords internal
digestionProductSummarization <- function(peptides, parameters) {
  message("#PEPTIDE SUMMARIZATION - Start\n")
  message(" + Summarization input:")
  message("  - A total number of ", nrow(peptides), " peptides is proceed for summarization.")

  # Create unique ID for each peptide based on the PTMType and PTMPos. No aggregation technique in any package supports lists...
  peptides$pep_id <- as.character(mapply(list, peptides$PTMType, peptides$PTMPos, SIMPLIFY = F))

  message("  - Unique peptide IDs are generated.")

  # Create a Sequence column where isoleucine is substituted by leucine.
  peptides$Sequence <- gsub("[I]", "L", peptides$Peptide)

  message("  - Isoleucine substitution to leucine is done.")

  # Helping functions for summarization per group for specific columns.
  log2.sum <- function(x) {
    x <- log2(sum(2^x, na.rm = T))
    if (is.finite(x)) {
      return(x)
    } else {
      return(NA)
    }
  }

  select.first <- function(x) {
    return(x[1])
  }

  # Create groups based on Sequence and pep_id
  peptides <- dplyr::group_by(.data = peptides, Sequence, pep_id)

  message("  - Peptide groups are generated.")

  peptides.1 <- peptides %>% dplyr::summarise(
    Peptide = list(Peptide),
    Start = list(Start),
    Stop = list(Stop),
    MC = list(MC),
    MZ1 = select.first(MZ1),
    MZ2 = select.first(MZ2),
    MZ3 = select.first(MZ3),
    Accession = list(Accession),
    Proteoform_ID = list(Proteoform_ID),
    PTMPos = select.first(PTMPos),
    PTMType = select.first(PTMType),
    Regulation_Amplitude = list(Regulation_Amplitude),
    Regulation_Pattern = list(Regulation_Pattern)
  )


  peptides.2 <- peptides %>% dplyr::summarise_at(.vars = dplyr::vars(parameters$QuantColnames), .funs = c("log2.sum"))

  peptides <- dplyr::inner_join(peptides.1, peptides.2, by = c("Sequence", "pep_id"))
  peptides <- dplyr::select(peptides, -c("pep_id"))

  message("  - Peptide groups summarization is done.")

  # Remove a percentage of randomly selected summarized peptides.
  remove <- sample(1:nrow(peptides), size = nrow(peptides) * parameters$LeastAbundantLoss, replace = FALSE)

  if (length(remove) != 0) {
    peptides <- peptides[-remove, ]
  }

  message("  - Remove ", parameters$LeastAbundantLoss * 100, "% of the least abundant peptides, which corresponds to ", length(remove), " peptides.")

  # add column with peptide detectability
  message("  - Calculating/predicting peptide detectability for later filtering with PeptideRanger.")
  peptides$Detectability <- addDetectability(peptides$Sequence, parameters)


  message(" + Summarization output:")
  message("  - A total number of ", nrow(peptides), " summarized peptides is generated.\n")
  message("#PEPTIDE SUMMARIZATION - Finish\n")

  return(peptides)
}
#####################

#' Create enriched and non-enriched fractions of proteolytic peptides
#'
#' This function separates modified peptides into enriched and non-enriched fractions,
#' adjusts the peptide abundances based on enrichment efficiency, and introduces noise
#' due to the enrichment process. It returns a list containing the enriched and
#' non-enriched peptide sets.
#'
#' @param DigestedProt A data frame containing the digested proteolytic peptides and associated data.
#' @param parameters A list of parameters that includes EnrichmentLoss, ModificationLoss, EnrichmentEfficiency,
#' EnrichmentNoise, and QuantColnames.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{NonEnriched}{A data frame containing the non-enriched peptide fraction, which includes both modified and non-modified peptides.}
#'   \item{Enriched}{A data frame containing the enriched peptide fraction, where modified peptides have been enriched based on the EnrichmentEfficiency, and noise has been added to simulate the enrichment process. If no modified peptides are present, this will be \code{NULL}.}
#' }
#'
filterDigestedProt <- function(DigestedProt, parameters) {
  modified <- lengths(DigestedProt$PTMType) != 0
  if (length(modified) == 0) modified <- NA

  ## Removing fraction according to ModificationLoss parameter
  numRemove <- 0
  if (sum(modified) > 0) {
    numRemove <- floor(sum(modified) * parameters$ModificationLoss)
  }
  message("\n#ENRICHMENT SIMULATION - Start\n")
  message(" + Modification loss")
  message("  - Remove ", numRemove, " modified peptides in non-enriched fraction according to parameter ModificationLoss (",
          parameters$ModificationLoss, ")")
  idx <- sample(which(modified), size = numRemove, replace = FALSE)
  nonenrichedtab <- DigestedProt
  if (length(idx) > 0)
    nonenrichedtab <- nonenrichedtab[-idx, ]

  if (sum(modified) == 0 | is.na(parameters$EnrichPTM) | parameters$EnrichmentEfficiency == 0) {
    message("\n#ENRICHMENT SIMULATION - Finish\n")
    return(list("NonEnriched" = nonenrichedtab, "Enriched" = NULL))
  } else {
    ## Exact copy of "sample"
    enrichedtab <- data.frame(DigestedProt)

    ## Removing fraction according to EnrichmentLoss parameter
    numRemove <- floor(nrow(enrichedtab) * parameters$EnrichmentLoss)
    message(" + Enrichment loss:")
    message("  - Remove ", numRemove, " peptides according to parameter EnrichmentLoss (", parameters$EnrichmentLoss, ")")
    idx <- sample(seq_len(nrow(enrichedtab)), size = numRemove, replace = FALSE)
    enrichedtab <- enrichedtab[-idx, ]

    # Select rows with PTM to be enriched
    modified <- sapply(enrichedtab$PTMType, function(x) sum(unlist(x) == parameters$EnrichPTM) > 0)

    message(" + Enriching PTM: ", parameters$EnrichPTM, ", having ", sum(modified), " peptides with this PTM in enriched fraction.")

    # Calculate total sum and average of intensities for modified and non-modified peptides in enriched fraction
    enrichedtab_modified <- enrichedtab[modified, parameters$QuantColnames]
    enrichedtab_nonmodified <- enrichedtab[!modified, parameters$QuantColnames]
    totalModified <- sum(unlist(enrichedtab_modified), na.rm = TRUE)
    averageModified <- mean(unlist(enrichedtab_modified), na.rm = TRUE)
    totalNonModified <- sum(unlist(enrichedtab_nonmodified), na.rm = TRUE)
    averageNonModified <- mean(unlist(enrichedtab_nonmodified), na.rm = TRUE)

    ## Adjust the intensities of modified peptides to mimic the mean difference of the parameter
    enrichedtab[!modified, parameters$QuantColnames] <- enrichedtab[!modified, parameters$QuantColnames] -
      parameters$EnrichmentQuantDiff
    enrichedtab[modified, parameters$QuantColnames] <- enrichedtab[modified, parameters$QuantColnames] +
      parameters$EnrichmentQuantDiff
    # Scale the total intensity of all peptides to the one before the adjustment
    enrichedtab[parameters$QuantColnames] <- enrichedtab[parameters$QuantColnames] *
      (totalModified + totalNonModified) / sum(unlist(enrichedtab[parameters$QuantColnames]), na.rm = TRUE)

    message(" + Enrichment efficiency:")
    message("  - Adjusted quantitative values to mimic the mean difference of the parameter EnrichmentQuantDiff (",
            parameters$EnrichmentQuantDiff, ").")
    # Calculate current fraction of modified peptides
    fracMod <- sum(modified) / nrow(enrichedtab)
    message("  - Modified peptides contribute to ",  fracMod * 100, "% of all present peptides in the enriched samples.")
    if (fracMod > parameters$EnrichmentEfficiency) {
      message("  - Warning: Enrichment efficiency is lower than the current fraction of modified peptides.
                No peptides will be removed")
    } else {
      # Remove the non-modified peptides from the enriched fraction to reach given enrichment efficiency
      numRemove <- floor(length(modified) - sum(modified) / parameters$EnrichmentEfficiency)
      message("  - Remove ", numRemove, " non-modified peptides to reach enrichment efficiency of ", parameters$EnrichmentEfficiency, "")
      idx <- sample(which(!modified), size = numRemove, replace = FALSE)
      enrichedtab <- enrichedtab[-idx, ]
      message("  - Enrichment efficiency is ", parameters$EnrichmentEfficiency, " leading to modified peptides contributing to ",
              parameters$EnrichmentEfficiency*100, "% of all present peptides in the enriched samples.")
    }

    ## Noise due to enrichment procedure:
    nrowTab <- nrow(enrichedtab)
    ncolTab <- length(parameters$QuantColnames)
    message(" + Enrichment noise:")
    message("  - The enrichment noise standard deviation is ", parameters$EnrichmentNoise, ".")
    mtx <- matrix(nrow = nrowTab, ncol = ncolTab, data = rnorm(n = ncolTab * nrowTab, mean = 0, sd = parameters$EnrichmentNoise))
    enrichedtab[, parameters$QuantColnames] <- enrichedtab[, parameters$QuantColnames] + mtx
    message("  - Noise added to all samples!")
    message("#ENRICHMENT SIMULATION - Finish")

    return(list("NonEnriched" = nonenrichedtab, "Enriched" = enrichedtab))
  }
}
#####################
