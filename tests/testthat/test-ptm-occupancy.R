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

# Build a minimal peptable with one modified and one unmodified row.
# mod_vals / unmod_vals are log2-scale vectors, one entry per column.
make_peptable <- function(mod_vals, unmod_vals,
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

  rbind(mod_row, unmod_row)
}

# ──────────────────────────────────────────────────────────────────────────────
# Basic correctness  (2 conditions × 2 replicates)
# ──────────────────────────────────────────────────────────────────────────────

test_that("calcPTMOccupancy returns one row per modified peptidoform", {
  pep <- make_peptable(mod_vals = c(1, 1, 2, 2), unmod_vals = c(1, 1, 1, 1))
  occ <- calcPTMOccupancy(pep, make_params())

  expect_equal(nrow(occ), 1L)
})

test_that("occupancy is 0.5 when modified and unmodified between-condition ratios are equal", {
  # Condition 1: mod=2, unmod=2;  Condition 2: mod=3, unmod=3
  # Rp = 2^(3-2)=2,  Ru = 2^(3-2)=2  =>  occ_C2 = 2/(2+2) = 0.5
  pep <- make_peptable(mod_vals = c(2, 2, 3, 3), unmod_vals = c(2, 2, 3, 3))
  occ <- calcPTMOccupancy(pep, make_params())

  expect_equal(as.numeric(occ[1, "C_2"]), 0.5, tolerance = 1e-9)
  expect_true(is.na(occ[1, "C_1"]))
})

test_that("occupancy approaches 1 when modified ratio >> unmodified ratio", {
  # Condition 1: mod=1, unmod=1;  Condition 2: mod=11, unmod=1
  # Rp = 2^10 = 1024,  Ru = 1  =>  occ_C2 = 1024/1025 ≈ 0.999
  pep <- make_peptable(mod_vals = c(1, 1, 11, 11), unmod_vals = c(1, 1, 1, 1))
  occ <- calcPTMOccupancy(pep, make_params())

  expect_true(as.numeric(occ[1, "C_2"]) > 0.99)
})

test_that("occupancy approaches 0 when modified ratio << unmodified ratio", {
  # Condition 1: mod=1, unmod=1;  Condition 2: mod=1, unmod=11
  # Rp = 1,  Ru = 2^10 = 1024  =>  occ_C2 = 1/1025 ≈ 0.001
  pep <- make_peptable(mod_vals = c(1, 1, 1, 1), unmod_vals = c(1, 1, 11, 11))
  occ <- calcPTMOccupancy(pep, make_params())

  expect_true(as.numeric(occ[1, "C_2"]) < 0.01)
})

# ──────────────────────────────────────────────────────────────────────────────
# Output structure
# ──────────────────────────────────────────────────────────────────────────────

test_that("output contains expected columns (Sequence, Accession, PTMPos, PTMType, C_1, C_2)", {
  pep <- make_peptable(mod_vals = c(1, 1, 2, 2), unmod_vals = c(1, 1, 1, 1))
  occ <- calcPTMOccupancy(pep, make_params())

  expect_true(all(c("Sequence", "Accession", "PTMPos", "PTMType", "C_1", "C_2") %in% names(occ)))
})

test_that("reference condition C_1 is always NA", {
  pep <- make_peptable(mod_vals = c(2, 2, 3, 3), unmod_vals = c(1, 1, 1, 1))
  occ <- calcPTMOccupancy(pep, make_params())

  expect_true(is.na(occ[1, "C_1"]))
})

test_that("non-reference occupancy values are in (0, 1)", {
  pep <- make_peptable(mod_vals = c(1, 1, 3, 3), unmod_vals = c(2, 2, 4, 4))
  occ <- calcPTMOccupancy(pep, make_params())

  val <- as.numeric(occ[1, "C_2"])
  expect_false(is.na(val))
  expect_true(val > 0 && val < 1)
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

test_that("NA in non-reference condition replicate propagates to that condition's occupancy only", {
  # 3 conditions, 1 replicate each: NA only in condition 2's modified replicate
  params <- make_params(num_cond = 3, num_reps = 1)
  qc     <- params$QuantColnames   # C_1_R_1, C_2_R_1, C_3_R_1

  pep <- make_peptable(
    mod_vals   = c(2, NA, 4),   # NA in condition 2 for mod
    unmod_vals = c(2, 2, 2),
    quant_cols = qc
  )
  occ <- calcPTMOccupancy(pep, params)

  expect_true(is.na(occ[1, "C_2"]))    # NA propagated to C_2
  expect_false(is.na(occ[1, "C_3"]))   # C_3 is unaffected
})

test_that("NA in reference condition replicate propagates to all non-reference conditions", {
  # NA in the reference condition's modified replicate → mod_mean_ref = NA
  # → log_Rp for all conditions = NA → all occupancies NA
  params <- make_params(num_cond = 2, num_reps = 2)
  pep <- make_peptable(
    mod_vals   = c(NA, 2, 3, 3),  # NA in C_1_R_1
    unmod_vals = c(2, 2, 2, 2)
  )
  occ <- calcPTMOccupancy(pep, params)

  expect_true(is.na(occ[1, "C_2"]))
})

test_that("NA in unmodified row in non-reference condition propagates to that condition", {
  params <- make_params(num_cond = 3, num_reps = 1)
  qc     <- params$QuantColnames

  mod_row <- data.frame(Sequence = "PEPTIDE", stringsAsFactors = FALSE)
  mod_row[qc] <- as.list(c(2, 3, 4))
  mod_row$PTMType   <- list("ph")
  mod_row$PTMPos    <- list(1L)
  mod_row$Accession <- list("P33333")

  # Two unmodified rows: second has NA in condition 2 only
  unmod1 <- data.frame(Sequence = "PEPTIDE", stringsAsFactors = FALSE)
  unmod1[qc] <- as.list(c(2, 2, 2))
  unmod1$PTMType <- list(character(0)); unmod1$PTMPos <- list(integer(0))
  unmod1$Accession <- list("P33333")

  unmod2 <- data.frame(Sequence = "PEPTIDE", stringsAsFactors = FALSE)
  unmod2[qc] <- as.list(c(2, NA, 2))   # NA in condition 2
  unmod2$PTMType <- list(character(0)); unmod2$PTMPos <- list(integer(0))
  unmod2$Accession <- list("P33333")

  pep <- rbind(mod_row, unmod1, unmod2)
  occ <- calcPTMOccupancy(pep, params)

  expect_true(is.na(occ[1, "C_2"]))    # C_2 affected
  expect_false(is.na(occ[1, "C_3"]))   # C_3 unaffected
})

# ──────────────────────────────────────────────────────────────────────────────
# Multiple unmodified rows: geometric mean per condition
# ──────────────────────────────────────────────────────────────────────────────

test_that("multiple unmodified rows are geometric-mean averaged per condition", {
  # 2 conditions, 1 replicate each
  params <- make_params(num_cond = 2, num_reps = 1)
  qc     <- params$QuantColnames   # C_1_R_1, C_2_R_1

  # Two unmod rows per condition: log2 values 0 and 2 in both conditions
  # Geometric mean = 2^((0+2)/2) = 2^1 = 2 in both conditions
  # Modified: 2 in both conditions (= geometric mean of unmod)
  # Rp = 2^(2-2) = 1,  Ru = 2^(2-2) = 1  =>  occ_C2 = 0.5
  mod_row <- data.frame(Sequence = "MYPEPTIDE", stringsAsFactors = FALSE)
  mod_row[qc] <- as.list(c(2, 2))
  mod_row$PTMType <- list("ph"); mod_row$PTMPos <- list(1L)
  mod_row$Accession <- list("P11111")

  unmod1 <- data.frame(Sequence = "MYPEPTIDE", stringsAsFactors = FALSE)
  unmod1[qc] <- as.list(c(0, 0))   # log2 = 0 in both conditions
  unmod1$PTMType <- list(character(0)); unmod1$PTMPos <- list(integer(0))
  unmod1$Accession <- list("P11111")

  unmod2 <- data.frame(Sequence = "MYPEPTIDE", stringsAsFactors = FALSE)
  unmod2[qc] <- as.list(c(2, 2))   # log2 = 2 in both conditions
  unmod2$PTMType <- list(character(0)); unmod2$PTMPos <- list(integer(0))
  unmod2$Accession <- list("P11111")

  pep <- rbind(mod_row, unmod1, unmod2)
  occ <- calcPTMOccupancy(pep, params)

  expect_equal(nrow(occ), 1L)
  expect_equal(as.numeric(occ[1, "C_2"]), 0.5, tolerance = 1e-9)
})

# ──────────────────────────────────────────────────────────────────────────────
# Multiple modified peptidoforms for the same sequence
# ──────────────────────────────────────────────────────────────────────────────

test_that("multiple modified rows for the same sequence each get their own row", {
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

  pep <- rbind(mod1, mod2, unmod)
  occ <- calcPTMOccupancy(pep, params)

  expect_equal(nrow(occ), 2L)
})

# ──────────────────────────────────────────────────────────────────────────────
# More than 2 conditions
# ──────────────────────────────────────────────────────────────────────────────

test_that("three conditions produce per-condition occupancy columns C_1, C_2, C_3", {
  params <- make_params(num_cond = 3, num_reps = 1)
  qc     <- params$QuantColnames   # C_1_R_1, C_2_R_1, C_3_R_1

  # Condition 1: mod=2, unmod=2 (reference)
  # Condition 2: mod=3, unmod=2  =>  Rp=2, Ru=1  =>  occ=2/3
  # Condition 3: mod=2, unmod=3  =>  Rp=1, Ru=2  =>  occ=1/3
  pep <- make_peptable(
    mod_vals   = c(2, 3, 2),
    unmod_vals = c(2, 2, 3),
    quant_cols = qc
  )
  occ <- calcPTMOccupancy(pep, params)

  expect_true(is.na(occ[1, "C_1"]))
  expect_equal(as.numeric(occ[1, "C_2"]), 2/3, tolerance = 1e-9)
  expect_equal(as.numeric(occ[1, "C_3"]), 1/3, tolerance = 1e-9)
})
