library(testthat)

# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────

# Build a parameters list with NumCond conditions × NumReps replicates.
# QuantColnames are ordered as: C_1_R_1, ..., C_1_R_n, C_2_R_1, ..., C_nc_R_n
make_params <- function(num_cond = 2, num_reps = 2) {
  quant_cols <- paste0(
    "C_", rep(seq_len(num_cond), each = num_reps),
    "_R_", rep(seq_len(num_reps), num_cond)
  )
  list(QuantColnames = quant_cols, NumCond = num_cond, NumReps = num_reps)
}

# Build a minimal peptable with one modified row, one counterpart unmodified row,
# and one non-counterpart unmodified row (different sequence, same protein).
# mod_vals / unmod_vals / other_vals are log2-scale vectors, one entry per column.
# other_vals defaults to unmod_vals for backward compatibility (Rprot = Ru when
# both share the same values, which is a degenerate case for the occupancy formula
# but sufficient for structural tests).
make_peptable <- function(mod_vals, unmod_vals, other_vals = unmod_vals,
                          ptmpos    = list(c(3L)),
                          ptmtype   = list(c("ph")),
                          quant_cols = make_params()$QuantColnames) {
  n_cols <- length(quant_cols)
  stopifnot(length(mod_vals) == n_cols, length(unmod_vals) == n_cols)

  mod_row <- data.frame(Sequence = "AASPEPR", stringsAsFactors = FALSE)
  mod_row[quant_cols] <- as.list(mod_vals)
  mod_row$PTMType    <- ptmtype
  mod_row$PTMPos     <- ptmpos
  mod_row$Accession  <- list("P12345")

  unmod_row <- data.frame(Sequence = "AASPEPR", stringsAsFactors = FALSE)
  unmod_row[quant_cols] <- as.list(unmod_vals)
  unmod_row$PTMType  <- list(character(0))
  unmod_row$PTMPos   <- list(integer(0))
  unmod_row$Accession <- list("P12345")

  # Non-counterpart unmodified peptide: different sequence, same protein.
  # Required for Rprot computation.
  other_row <- data.frame(Sequence = "BGPEPTIDE", stringsAsFactors = FALSE)
  other_row[quant_cols] <- as.list(other_vals)
  other_row$PTMType  <- list(character(0))
  other_row$PTMPos   <- list(integer(0))
  other_row$Accession <- list("P12345")

  rbind(mod_row, unmod_row, other_row)
}

# ──────────────────────────────────────────────────────────────────────────────
# Basic correctness  (2 conditions × 2 replicates)
# ──────────────────────────────────────────────────────────────────────────────

test_that("calcPTMOccupancy returns one row per modified peptidoform", {
  pep <- make_peptable(mod_vals = c(1, 1, 2, 2), unmod_vals = c(1, 1, 1, 1))
  occ <- calcPTMOccupancy(pep, make_params())

  expect_equal(nrow(occ), 1L)
})

test_that("calcPTMOccupancy skips peptides with ambiguous protein start positions", {
  pep <- make_peptable(mod_vals = c(1, 1, 2, 2), unmod_vals = c(1, 1, 1, 1))
  pep$Start <- I(list(c(10L, 20L), 10L, 40L))

  occ <- calcPTMOccupancy(pep, make_params())

  expect_equal(nrow(occ), 0L)
})

test_that("calcPTMOccupancy skips peptides with ambiguous accessions", {
  pep <- make_peptable(mod_vals = c(1, 1, 2, 2), unmod_vals = c(1, 1, 1, 1))
  pep$Accession[[1]] <- c("P12345", "P67890")

  occ <- calcPTMOccupancy(pep, make_params())

  expect_equal(nrow(occ), 0L)
})

test_that("formula correctly estimates reference-condition occupancy occ_1", {
  # 2 conditions, 1 replicate each.
  # mod: ref=0, C2=2 (log2)  =>  Rm = 2^(2-0) = 4
  # unmod: ref=0, C2=0        =>  Ru = 2^(0-0) = 1
  # other: ref=0, C2=1        =>  Rprot = 2^(1-0) = 2
  # occ_1 = (Rprot - Ru) / (Rm - Ru) = (2 - 1) / (4 - 1) = 1/3
  params <- make_params(num_cond = 2, num_reps = 1)
  qc     <- params$QuantColnames
  pep    <- make_peptable(mod_vals   = c(0, 2), unmod_vals = c(0, 0),
                          other_vals = c(0, 1), quant_cols = qc)
  occ    <- calcPTMOccupancy(pep, params)

  expect_equal(as.numeric(occ[1, "C_1"]), 1/3, tolerance = 1e-9)
})

test_that("formula correctly computes per-condition occupancy as occ_1 * Rm / Rprot", {
  # Same setup as above: occ_1 = 1/3, Rm = 4, Rprot = 2
  # occ_2 = occ_1 * Rm / Rprot = (1/3) * 4 / 2 = 2/3
  params <- make_params(num_cond = 2, num_reps = 1)
  qc     <- params$QuantColnames
  pep    <- make_peptable(mod_vals   = c(0, 2), unmod_vals = c(0, 0),
                          other_vals = c(0, 1), quant_cols = qc)
  occ    <- calcPTMOccupancy(pep, params)

  expect_equal(as.numeric(occ[1, "C_2"]), 2/3, tolerance = 1e-9)
})

test_that("occupancy approaches 0 when modified ratio << protein ratio", {
  # 2 conditions, 2 replicates: mod unchanged (Rm=1), unmod and protein increase 1024x (Ru=Rprot=1024).
  # occ_1 = (Rprot - Ru) / (Rm - Ru) = (1024-1024)/(1-1024) = 0
  # occ_2 = occ_1 * Rm / Rprot = 0  (denominator is degenerate: Rm=Ru)
  pep <- make_peptable(mod_vals = c(1, 1, 1, 1), unmod_vals = c(1, 1, 11, 11))
  occ <- calcPTMOccupancy(pep, make_params())

  expect_true(as.numeric(occ[1, "C_2"]) < 0.01)
})

# ──────────────────────────────────────────────────────────────────────────────
# Output structure
# ──────────────────────────────────────────────────────────────────────────────

test_that("output contains expected columns (Sequence, Accession, PTMPos, ProteinPTMPos, PTMType, C_1, C_2)", {
  pep <- make_peptable(mod_vals = c(1, 1, 2, 2), unmod_vals = c(1, 1, 1, 1))
  occ <- calcPTMOccupancy(pep, make_params())

  expect_true(all(c("Sequence", "Accession", "PTMPos", "ProteinPTMPos", "PTMType", "C_1", "C_2") %in% names(occ)))
})

test_that("ProteinPTMPos reports protein-level modification positions when Start is available", {
  pep <- make_peptable(mod_vals = c(1, 1, 2, 2), unmod_vals = c(1, 1, 1, 1))
  pep$Start <- I(list(10L, 10L, 40L))

  occ <- calcPTMOccupancy(pep, make_params())

  expect_identical(occ$PTMPos[[1]], 3L)
  expect_identical(occ$ProteinPTMPos[[1]], 12L)
})

test_that("C_1 column contains the estimated reference-condition occupancy occ_1", {
  # Decreasing occupancy case: mod stays, unmod increases, protein doubles.
  # mod: ref=0, C2=0  =>  Rm = 1
  # unmod: ref=0, C2=2 =>  Ru = 4
  # other: ref=0, C2=1 =>  Rprot = 2
  # occ_1 = (Rprot - Ru) / (Rm - Ru) = (2 - 4) / (1 - 4) = 2/3
  params <- make_params(num_cond = 2, num_reps = 1)
  qc     <- params$QuantColnames
  pep    <- make_peptable(mod_vals   = c(0, 0), unmod_vals = c(0, 2),
                          other_vals = c(0, 1), quant_cols = qc)
  occ    <- calcPTMOccupancy(pep, params)

  expect_equal(as.numeric(occ[1, "C_1"]), 2/3, tolerance = 1e-9)
})

test_that("both C_1 and non-reference occupancy values are in (0, 1) for consistent data", {
  # Canonical consistent data: occ_1 = 1/3, occ_2 = 2/3
  params <- make_params(num_cond = 2, num_reps = 1)
  qc     <- params$QuantColnames
  pep    <- make_peptable(mod_vals   = c(0, 2), unmod_vals = c(0, 0),
                          other_vals = c(0, 1), quant_cols = qc)
  occ    <- calcPTMOccupancy(pep, params)

  c1_val <- as.numeric(occ[1, "C_1"])
  c2_val <- as.numeric(occ[1, "C_2"])
  expect_false(is.na(c1_val))
  expect_true(c1_val > 0 && c1_val < 1)
  expect_false(is.na(c2_val))
  expect_true(c2_val > 0 && c2_val < 1)
})

# ──────────────────────────────────────────────────────────────────────────────
# Edge cases
# ──────────────────────────────────────────────────────────────────────────────

test_that("returns empty data.frame when no modified peptides exist", {
  params <- make_params()
  qc     <- params$QuantColnames

  unmod_row <- data.frame(Sequence = "PEPTIDE", stringsAsFactors = FALSE)
  unmod_row[qc] <- as.list(rep(2, length(qc)))
  unmod_row$PTMType   <- list(character(0))
  unmod_row$PTMPos    <- list(integer(0))
  unmod_row$Accession <- list("P99999")

  occ <- calcPTMOccupancy(unmod_row, params)
  expect_equal(nrow(occ), 0L)
})

test_that("returns empty data.frame when modified peptide has no unmodified counterpart", {
  params <- make_params()
  qc     <- params$QuantColnames

  mod_row <- data.frame(Sequence = "UNIQUEPEP", stringsAsFactors = FALSE)
  mod_row[qc] <- as.list(rep(3, length(qc)))
  mod_row$PTMType   <- list("ph")
  mod_row$PTMPos    <- list(2L)
  mod_row$Accession <- list("P00000")

  occ <- calcPTMOccupancy(mod_row, params)
  expect_equal(nrow(occ), 0L)
})

test_that("returns empty data.frame when NumCond < 2", {
  pep <- make_peptable(
    mod_vals   = c(1, 2),
    unmod_vals = c(1, 2),
    quant_cols = c("C_1_R_1", "C_1_R_2")
  )
  occ <- calcPTMOccupancy(pep, list(QuantColnames = c("C_1_R_1", "C_1_R_2"),
                                    NumCond = 1L, NumReps = 2L))
  expect_equal(nrow(occ), 0L)
})

test_that("returns empty data.frame when NumCond/NumReps absent from parameters", {
  pep <- make_peptable(mod_vals = c(1, 1, 2, 2), unmod_vals = c(1, 1, 1, 1))
  occ <- calcPTMOccupancy(pep, list(QuantColnames = make_params()$QuantColnames))
  expect_equal(nrow(occ), 0L)
})

# ──────────────────────────────────────────────────────────────────────────────
# NA propagation
# ──────────────────────────────────────────────────────────────────────────────

test_that("NA in a single-replicate condition propagates via occ_1 to all output columns", {
  # occ_1 is the mean of per-condition terms.  A NaN in any condition (arising
  # from mean(c(NA), na.rm = TRUE) = NaN) propagates through the mean because
  # NaN is not NA and is not removed by na.rm = TRUE.  The result is occ_1 = NaN,
  # which makes all columns NaN.  In R, is.na(NaN) == TRUE.
  params <- make_params(num_cond = 3, num_reps = 1)
  qc     <- params$QuantColnames   # C_1_R_1, C_2_R_1, C_3_R_1

  pep <- make_peptable(
    mod_vals   = c(0, NA, 2),   # NA in condition 2's only replicate → mod_mean[2] = NaN
    unmod_vals = c(0,  0, 0),
    other_vals = c(0,  1, 1),
    quant_cols = qc
  )
  occ <- calcPTMOccupancy(pep, params)

  expect_true(is.na(occ[1, "C_1"]))  # occ_1 = NaN (is.na(NaN) == TRUE)
  expect_true(is.na(occ[1, "C_2"]))  # NaN propagated via occ_1
  expect_true(is.na(occ[1, "C_3"]))  # NaN propagated via occ_1
})

test_that("NA in one of multiple replicates is silently removed (na.rm = TRUE)", {
  # mean(c(NA, x), na.rm = TRUE) = x: a single missing replicate is excluded and
  # the remaining replicate(s) are used.  The result equals the no-NA case.
  params <- make_params(num_cond = 2, num_reps = 2)
  qc     <- params$QuantColnames   # C_1_R_1, C_1_R_2, C_2_R_1, C_2_R_2

  pep_with_na <- make_peptable(
    mod_vals   = c(NA, 0, 2, 2),  # NA in C_1_R_1; mean of C_1 = 0
    unmod_vals = c(0, 0, 0, 0),
    other_vals = c(0, 0, 1, 1),
    quant_cols = qc
  )
  pep_complete <- make_peptable(
    mod_vals   = c(0, 0, 2, 2),   # same values without NA
    unmod_vals = c(0, 0, 0, 0),
    other_vals = c(0, 0, 1, 1),
    quant_cols = qc
  )

  occ_na  <- calcPTMOccupancy(pep_with_na,  params)
  occ_cmp <- calcPTMOccupancy(pep_complete, params)

  expect_equal(as.numeric(occ_na[1, "C_1"]), as.numeric(occ_cmp[1, "C_1"]), tolerance = 1e-9)
  expect_equal(as.numeric(occ_na[1, "C_2"]), as.numeric(occ_cmp[1, "C_2"]), tolerance = 1e-9)
})

test_that("NA in non-counterpart background propagates via occ_1 to all output columns", {
  # Rprot for the NA condition becomes NaN, which propagates through the occ_1
  # mean to all columns (NaN is not removed by na.rm = TRUE in mean()).
  params <- make_params(num_cond = 3, num_reps = 1)
  qc     <- params$QuantColnames

  mod_row <- data.frame(Sequence = "PEPTIDE", stringsAsFactors = FALSE)
  mod_row[qc] <- as.list(c(0, 2, 2))
  mod_row$PTMType   <- list("ph")
  mod_row$PTMPos    <- list(1L)
  mod_row$Accession <- list("P33333")

  # Counterpart unmodified (same sequence as modified)
  unmod_same <- data.frame(Sequence = "PEPTIDE", stringsAsFactors = FALSE)
  unmod_same[qc] <- as.list(c(0, 0, 0))
  unmod_same$PTMType <- list(character(0)); unmod_same$PTMPos <- list(integer(0))
  unmod_same$Accession <- list("P33333")

  # Non-counterpart unmodified: NA in condition 2 → Rprot_2 = NaN
  other_row <- data.frame(Sequence = "OTHERPEP", stringsAsFactors = FALSE)
  other_row[qc] <- as.list(c(0, NA, 1))
  other_row$PTMType <- list(character(0)); other_row$PTMPos <- list(integer(0))
  other_row$Accession <- list("P33333")

  pep <- rbind(mod_row, unmod_same, other_row)
  occ <- calcPTMOccupancy(pep, params)

  expect_true(is.na(occ[1, "C_1"]))   # occ_1 = NaN (is.na(NaN) == TRUE)
  expect_true(is.na(occ[1, "C_2"]))   # NaN propagated via occ_1
  expect_true(is.na(occ[1, "C_3"]))   # NaN propagated via occ_1
})

# ──────────────────────────────────────────────────────────────────────────────
# Multiple unmodified rows: geometric mean per condition
# ──────────────────────────────────────────────────────────────────────────────

test_that("sequence with multiple unmodified counterpart rows is skipped with a warning", {
  # The code requires exactly one unmodified row per sequence.  When two unmodified
  # rows exist for the same sequence a warning is issued and the sequence is skipped.
  params <- make_params(num_cond = 2, num_reps = 1)
  qc     <- params$QuantColnames   # C_1_R_1, C_2_R_1

  mod_row <- data.frame(Sequence = "MYPEPTIDE", stringsAsFactors = FALSE)
  mod_row[qc] <- as.list(c(2, 2))
  mod_row$PTMType <- list("ph"); mod_row$PTMPos <- list(1L)
  mod_row$Accession <- list("P11111")

  unmod1 <- data.frame(Sequence = "MYPEPTIDE", stringsAsFactors = FALSE)
  unmod1[qc] <- as.list(c(1, 1))
  unmod1$PTMType <- list(character(0)); unmod1$PTMPos <- list(integer(0))
  unmod1$Accession <- list("P11111")

  unmod2 <- data.frame(Sequence = "MYPEPTIDE", stringsAsFactors = FALSE)
  unmod2[qc] <- as.list(c(3, 3))
  unmod2$PTMType <- list(character(0)); unmod2$PTMPos <- list(integer(0))
  unmod2$Accession <- list("P11111")

  # Non-counterpart unmodified row required for Rprot
  other <- data.frame(Sequence = "OTHERPEP", stringsAsFactors = FALSE)
  other[qc] <- as.list(c(2, 3))
  other$PTMType <- list(character(0)); other$PTMPos <- list(integer(0))
  other$Accession <- list("P11111")

  pep <- rbind(mod_row, unmod1, unmod2, other)

  expect_warning(
    occ <- calcPTMOccupancy(pep, params),
    regexp = "more than one peptide"
  )
  expect_equal(nrow(occ), 0L)
})

# ──────────────────────────────────────────────────────────────────────────────
# Multiple modified peptidoforms for the same sequence
# ──────────────────────────────────────────────────────────────────────────────

test_that("sequence with multiple modified peptidoforms is skipped with a warning", {
  # The code requires exactly one modified row per sequence.  When two modified
  # rows exist for the same sequence a warning is issued and the sequence is skipped.
  params <- make_params(num_cond = 2, num_reps = 1)
  qc     <- params$QuantColnames

  mod1 <- data.frame(Sequence = "STYPEP", stringsAsFactors = FALSE)
  mod1[qc] <- as.list(c(2, 3))
  mod1$PTMType <- list("ph"); mod1$PTMPos <- list(1L)
  mod1$Accession <- list("P22222")

  mod2 <- data.frame(Sequence = "STYPEP", stringsAsFactors = FALSE)
  mod2[qc] <- as.list(c(2, 4))
  mod2$PTMType <- list("ph"); mod2$PTMPos <- list(2L)
  mod2$Accession <- list("P22222")

  unmod <- data.frame(Sequence = "STYPEP", stringsAsFactors = FALSE)
  unmod[qc] <- as.list(c(2, 2))
  unmod$PTMType <- list(character(0)); unmod$PTMPos <- list(integer(0))
  unmod$Accession <- list("P22222")

  # Non-counterpart unmodified row required for Rprot
  other <- data.frame(Sequence = "OTHERPEP", stringsAsFactors = FALSE)
  other[qc] <- as.list(c(2, 3))
  other$PTMType <- list(character(0)); other$PTMPos <- list(integer(0))
  other$Accession <- list("P22222")

  pep <- rbind(mod1, mod2, unmod, other)

  expect_warning(
    occ <- calcPTMOccupancy(pep, params),
    regexp = "more than one peptide"
  )
  expect_equal(nrow(occ), 0L)
})

# ──────────────────────────────────────────────────────────────────────────────
# More than 2 conditions
# ──────────────────────────────────────────────────────────────────────────────

test_that("three conditions produce per-condition occupancy columns C_1, C_2, C_3", {
  # 3 conditions, 1 replicate each.
  # mod: ref=0, C2=2, C3=2  =>  Rm = [4, 4]
  # unmod: ref=0, C2=0, C3=0 =>  Ru = [1, 1]
  # other: ref=0, C2=1, C3=1 =>  Rprot = [2, 2]
  # occ_1 = mean([(2-1)/(4-1), (2-1)/(4-1)]) = 1/3
  # occ_2 = (1/3) * 4 / 2 = 2/3
  # occ_3 = (1/3) * 4 / 2 = 2/3
  params <- make_params(num_cond = 3, num_reps = 1)
  qc     <- params$QuantColnames
  pep    <- make_peptable(mod_vals   = c(0, 2, 2),
                          unmod_vals = c(0, 0, 0),
                          other_vals = c(0, 1, 1),
                          quant_cols = qc)
  occ    <- calcPTMOccupancy(pep, params)

  expect_equal(as.numeric(occ[1, "C_1"]), 1/3, tolerance = 1e-9)
  expect_equal(as.numeric(occ[1, "C_2"]), 2/3, tolerance = 1e-9)
  expect_equal(as.numeric(occ[1, "C_3"]), 2/3, tolerance = 1e-9)
})

test_that("formula correctly handles decreasing occupancy (occ_1 > occ_2)", {
  # Decreasing occupancy case where occ_1 = 2/3 and occ_2 = 1/3.
  # mod: ref=0, C2=0  =>  Rm = 2^(0-0) = 1
  # unmod: ref=0, C2=2 =>  Ru = 2^(2-0) = 4
  # other: ref=0, C2=1 =>  Rprot = 2^(1-0) = 2
  # occ_1 = (2 - 4) / (1 - 4) = (-2)/(-3) = 2/3
  # occ_2 = (2/3) * 1 / 2 = 1/3
  params <- make_params(num_cond = 2, num_reps = 1)
  qc     <- params$QuantColnames

  pep <- make_peptable(mod_vals   = c(0, 0),
                       unmod_vals = c(0, 2),
                       other_vals = c(0, 1),
                       quant_cols = qc)
  occ <- calcPTMOccupancy(pep, params)

  expect_equal(as.numeric(occ[1, "C_1"]), 2/3, tolerance = 1e-9)
  expect_equal(as.numeric(occ[1, "C_2"]), 1/3, tolerance = 1e-9)
})

# ──────────────────────────────────────────────────────────────────────────────
# Protein ratio from non-counterpart unmodified peptides only (requirement c)
# ──────────────────────────────────────────────────────────────────────────────

test_that("protein ratio uses only non-counterpart unmodified peptides", {
  # Protein P12345 has two peptides:
  #   - MODSEQ: observed as both modified and unmodified (counterpart)
  #   - OTHERPEP: observed as unmodified only (non-counterpart)
  # The protein ratio Rprot must use ONLY OTHERPEP; the counterpart is excluded.
  #
  # 2 conditions, 1 replicate each.
  # Modified MODSEQ:       ref=0, C2=2  =>  Rm = 2^(2-0) = 4
  # Unmodified MODSEQ:     ref=0, C2=0  =>  Ru = 1  (excluded from Rprot)
  # Unmodified OTHERPEP:   ref=0, C2=1  =>  Rprot = 2^(1-0) = 2
  # occ_1 = (2 - 1) / (4 - 1) = 1/3
  # occ_2 = (1/3) * 4 / 2 = 2/3
  # If counterpart (ref=0,C2=0) were also included in Rprot:
  #   colMeans of [0,0] and [0,1] = [0, 0.5], Rprot = sqrt(2) ≈ 1.414
  #   occ_1 ≈ (1.414 - 1)/(4 - 1) ≈ 0.138, occ_2 ≈ 0.39  (different!)
  params <- make_params(num_cond = 2, num_reps = 1)
  qc     <- params$QuantColnames

  mod_row <- data.frame(Sequence = "MODSEQ", stringsAsFactors = FALSE)
  mod_row[qc] <- as.list(c(0, 2))
  mod_row$PTMType   <- list("ph")
  mod_row$PTMPos    <- list(1L)
  mod_row$Accession <- list("P12345")

  unmod_same <- data.frame(Sequence = "MODSEQ", stringsAsFactors = FALSE)
  unmod_same[qc] <- as.list(c(0, 0))
  unmod_same$PTMType   <- list(character(0))
  unmod_same$PTMPos    <- list(integer(0))
  unmod_same$Accession <- list("P12345")

  unmod_other <- data.frame(Sequence = "OTHERPEP", stringsAsFactors = FALSE)
  unmod_other[qc] <- as.list(c(0, 1))
  unmod_other$PTMType   <- list(character(0))
  unmod_other$PTMPos    <- list(integer(0))
  unmod_other$Accession <- list("P12345")

  pep <- rbind(mod_row, unmod_same, unmod_other)
  occ <- calcPTMOccupancy(pep, params)

  expect_equal(nrow(occ), 1L)
  # Must equal 2/3 (Rprot from OTHERPEP only), not ≈ 0.39 (if counterpart included)
  expect_equal(as.numeric(occ[1, "C_2"]), 2/3, tolerance = 1e-9)
})

test_that("returns empty data.frame when modified peptide has no non-counterpart unmodified peptides", {
  # Only a counterpart unmodified exists (same sequence); no other unmodified
  # peptides from the same protein → Rprot cannot be computed → skip.
  params <- make_params()
  qc     <- params$QuantColnames

  mod_row <- data.frame(Sequence = "UNIQUEPEP", stringsAsFactors = FALSE)
  mod_row[qc] <- as.list(rep(3, length(qc)))
  mod_row$PTMType   <- list("ph")
  mod_row$PTMPos    <- list(2L)
  mod_row$Accession <- list("P00001")

  unmod_row <- data.frame(Sequence = "UNIQUEPEP", stringsAsFactors = FALSE)
  unmod_row[qc] <- as.list(rep(2, length(qc)))
  unmod_row$PTMType  <- list(character(0))
  unmod_row$PTMPos   <- list(integer(0))
  unmod_row$Accession <- list("P00001")

  occ <- calcPTMOccupancy(rbind(mod_row, unmod_row), params)
  expect_equal(nrow(occ), 0L)
})
