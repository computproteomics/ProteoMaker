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
# mod_vals / unmod_vals are log2-scale vectors, one entry per column.
# The non-counterpart row uses the same values as unmod_vals so that protein-
# ratio (Rprot) behaviour matches unmod_mean in tests that do not specifically
# exercise the two values being different.
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

  # Non-counterpart unmodified peptide: different sequence, same protein.
  # Required for Rprot computation (requirement c).
  other_row <- data.frame(Sequence = "BGPEPTIDE", stringsAsFactors = FALSE)
  other_row[quant_cols] <- as.list(unmod_vals)
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

test_that("occupancy equals reference occupancy when modified and unmodified between-condition ratios are equal", {
  # Condition 1: mod=2, unmod=2  =>  occ_ref = 4/(4+4) = 0.5
  # Condition 2: mod=3, unmod=3  =>  Rp=2, Ru=2
  # When Rp == Ru, occ_c = occ_ref  =>  occ_C2 = 0.5
  pep <- make_peptable(mod_vals = c(2, 2, 3, 3), unmod_vals = c(2, 2, 3, 3))
  occ <- calcPTMOccupancy(pep, make_params())

  expect_equal(as.numeric(occ[1, "C_2"]), 0.5, tolerance = 1e-9)
  expect_true(is.na(occ[1, "C_1"]))
})

test_that("occupancy approaches 1 when modified ratio >> protein ratio", {
  # Condition 1: mod=1, unmod=1  =>  occ_ref = 0.5
  # Condition 2: mod=11, unmod=1  =>  Rp=2^10=1024, Rprot=1 (one unmod peptide)
  # occ_C2 = (0.5*1024)/(0.5*1024 + 0.5*1) = 512/512.5 ≈ 0.999
  pep <- make_peptable(mod_vals = c(1, 1, 11, 11), unmod_vals = c(1, 1, 1, 1))
  occ <- calcPTMOccupancy(pep, make_params())

  expect_true(as.numeric(occ[1, "C_2"]) > 0.99)
})

test_that("occupancy approaches 0 when modified ratio << protein ratio", {
  # Condition 1: mod=1, unmod=1  =>  occ_ref = 0.5
  # Condition 2: mod=1, unmod=11  =>  Rp=1, Rprot=2^10=1024 (one unmod peptide)
  # occ_C2 = (0.5*1)/(0.5*1 + 0.5*1024) = 0.5/512.5 ≈ 0.001
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

test_that("NA in non-counterpart unmodified row in non-reference condition propagates to that condition", {
  params <- make_params(num_cond = 3, num_reps = 1)
  qc     <- params$QuantColnames

  mod_row <- data.frame(Sequence = "PEPTIDE", stringsAsFactors = FALSE)
  mod_row[qc] <- as.list(c(2, 3, 4))
  mod_row$PTMType   <- list("ph")
  mod_row$PTMPos    <- list(1L)
  mod_row$Accession <- list("P33333")

  # Counterpart unmodified (same sequence as modified)
  unmod_same <- data.frame(Sequence = "PEPTIDE", stringsAsFactors = FALSE)
  unmod_same[qc] <- as.list(c(2, 2, 2))
  unmod_same$PTMType <- list(character(0)); unmod_same$PTMPos <- list(integer(0))
  unmod_same$Accession <- list("P33333")

  # Non-counterpart unmodified: NA in condition 2 only → affects Rprot for C_2
  other_row <- data.frame(Sequence = "OTHERPEP", stringsAsFactors = FALSE)
  other_row[qc] <- as.list(c(2, NA, 2))
  other_row$PTMType <- list(character(0)); other_row$PTMPos <- list(integer(0))
  other_row$Accession <- list("P33333")

  pep <- rbind(mod_row, unmod_same, other_row)
  occ <- calcPTMOccupancy(pep, params)

  expect_true(is.na(occ[1, "C_2"]))    # C_2 affected by NA in non-counterpart
  expect_false(is.na(occ[1, "C_3"]))   # C_3 unaffected
})

# ──────────────────────────────────────────────────────────────────────────────
# Multiple unmodified rows: geometric mean per condition
# ──────────────────────────────────────────────────────────────────────────────

test_that("multiple counterpart unmodified rows are geometric-mean averaged for occ_ref", {
  # 2 conditions, 1 replicate each
  params <- make_params(num_cond = 2, num_reps = 1)
  qc     <- params$QuantColnames   # C_1_R_1, C_2_R_1

  # Two counterpart unmod rows: log2 = 1 and 3 in both conditions.
  # colMeans in log2 space = (1+3)/2 = 2  <=>  geometric mean in linear = 2^2 = 4.
  # Arithmetic mean of linear values would give (2+8)/2 = 5, i.e. log2 ≈ 2.32 — different!
  # Modified: log2 = 2 in both conditions (equals the geometric-mean unmod).
  # occ_ref = 2^2 / (2^2 + 2^2) = 0.5
  # Rp = 2^(2-2) = 1, Rprot (from OTHER, same values) = 1  =>  occ_C2 = 0.5
  mod_row <- data.frame(Sequence = "MYPEPTIDE", stringsAsFactors = FALSE)
  mod_row[qc] <- as.list(c(2, 2))
  mod_row$PTMType <- list("ph"); mod_row$PTMPos <- list(1L)
  mod_row$Accession <- list("P11111")

  unmod1 <- data.frame(Sequence = "MYPEPTIDE", stringsAsFactors = FALSE)
  unmod1[qc] <- as.list(c(1, 1))   # log2 = 1 in both conditions
  unmod1$PTMType <- list(character(0)); unmod1$PTMPos <- list(integer(0))
  unmod1$Accession <- list("P11111")

  unmod2 <- data.frame(Sequence = "MYPEPTIDE", stringsAsFactors = FALSE)
  unmod2[qc] <- as.list(c(3, 3))   # log2 = 3 in both conditions
  unmod2$PTMType <- list(character(0)); unmod2$PTMPos <- list(integer(0))
  unmod2$Accession <- list("P11111")

  # Non-counterpart unmodified row required for Rprot (same values → Rprot = 1)
  other <- data.frame(Sequence = "OTHERPEP", stringsAsFactors = FALSE)
  other[qc] <- as.list(c(2, 2))
  other$PTMType <- list(character(0)); other$PTMPos <- list(integer(0))
  other$Accession <- list("P11111")

  pep <- rbind(mod_row, unmod1, unmod2, other)
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

  # Non-counterpart unmodified row required for Rprot
  other <- data.frame(Sequence = "OTHERPEP", stringsAsFactors = FALSE)
  other[qc] <- as.list(c(2, 2))
  other$PTMType <- list(character(0)); other$PTMPos <- list(integer(0))
  other$Accession <- list("P22222")

  pep <- rbind(mod1, mod2, unmod, other)
  occ <- calcPTMOccupancy(pep, params)

  expect_equal(nrow(occ), 2L)
})

# ──────────────────────────────────────────────────────────────────────────────
# More than 2 conditions
# ──────────────────────────────────────────────────────────────────────────────

test_that("three conditions produce per-condition occupancy columns C_1, C_2, C_3", {
  params <- make_params(num_cond = 3, num_reps = 1)
  qc     <- params$QuantColnames   # C_1_R_1, C_2_R_1, C_3_R_1

  # Condition 1: mod=2, unmod=2  =>  occ_ref = 4/(4+4) = 0.5
  # Condition 2: mod=3, unmod=2  =>  Rp=2, Rprot=1  =>  occ=(0.5*2)/(0.5*2+0.5*1)=2/3
  # Condition 3: mod=2, unmod=3  =>  Rp=1, Rprot=2  =>  occ=(0.5*1)/(0.5*1+0.5*2)=1/3
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

test_that("3-ratio formula uses reference occupancy: result differs from 2-ratio when occ_ref != 0.5", {
  # When occ_ref = 0.5 the formulas coincide; here we use occ_ref != 0.5 to
  # verify the full formula is applied.
  # Condition 1: mod=4 (linear 16), unmod=2 (linear 4)
  #   occ_ref = 16/(16+4) = 0.8
  # Condition 2: mod=4 (same), unmod=4 (linear 16)
  #   Rp = 2^(4-4) = 1,  Rprot = 2^(4-2) = 4
  # 3-ratio: occ_C2 = (0.8*1)/(0.8*1+0.2*4) = 0.8/1.6 = 0.5
  # 2-ratio (wrong): would give 1/(1+4) = 0.2
  params <- make_params(num_cond = 2, num_reps = 1)
  qc     <- params$QuantColnames

  pep <- make_peptable(
    mod_vals   = c(4, 4),
    unmod_vals = c(2, 4),
    quant_cols = qc
  )
  occ <- calcPTMOccupancy(pep, params)

  expect_equal(as.numeric(occ[1, "C_2"]), 0.5, tolerance = 1e-9)
})

# ──────────────────────────────────────────────────────────────────────────────
# Protein ratio from non-counterpart unmodified peptides only (requirement c)
# ──────────────────────────────────────────────────────────────────────────────

test_that("protein ratio uses only non-counterpart unmodified peptides", {
  # Protein P12345 has two peptides:
  #   - MODSEQ: observed as both modified and unmodified (counterpart)
  #   - OTHERPEP: only observed as unmodified (non-counterpart)
  # The protein ratio Rprot should use ONLY OTHERPEP, excluding MODSEQ unmodified.
  #
  # 2 conditions, 1 replicate each.
  # Modified MODSEQ:       C1=2, C2=3  =>  Rp = 2^(3-2) = 2
  # Unmodified MODSEQ:     C1=2, C2=2  =>  counterpart (excluded from Rprot)
  # Unmodified OTHERPEP:   C1=2, C2=4  =>  non-counterpart: prot_mean = c(2,4), Rprot = 4
  # occ_ref = 2^2/(2^2+2^2) = 0.5
  # occ_C2 = (0.5*2)/(0.5*2+0.5*4) = 1/(1+2) = 1/3
  # If counterpart also included (wrong): prot_mean[2]=(2+4)/2=3, Rprot=2 => occ_C2=0.5
  params <- make_params(num_cond = 2, num_reps = 1)
  qc     <- params$QuantColnames

  mod_row <- data.frame(Sequence = "MODSEQ", stringsAsFactors = FALSE)
  mod_row[qc] <- as.list(c(2, 3))
  mod_row$PTMType   <- list("ph")
  mod_row$PTMPos    <- list(1L)
  mod_row$Accession <- list("P12345")

  unmod_same <- data.frame(Sequence = "MODSEQ", stringsAsFactors = FALSE)
  unmod_same[qc] <- as.list(c(2, 2))
  unmod_same$PTMType   <- list(character(0))
  unmod_same$PTMPos    <- list(integer(0))
  unmod_same$Accession <- list("P12345")

  unmod_other <- data.frame(Sequence = "OTHERPEP", stringsAsFactors = FALSE)
  unmod_other[qc] <- as.list(c(2, 4))
  unmod_other$PTMType   <- list(character(0))
  unmod_other$PTMPos    <- list(integer(0))
  unmod_other$Accession <- list("P12345")

  pep <- rbind(mod_row, unmod_same, unmod_other)
  occ <- calcPTMOccupancy(pep, params)

  expect_equal(nrow(occ), 1L)
  # Must equal 1/3 (Rprot from OTHERPEP only), NOT 0.5 (which would include counterpart)
  expect_equal(as.numeric(occ[1, "C_2"]), 1/3, tolerance = 1e-9)
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
