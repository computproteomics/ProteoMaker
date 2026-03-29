library(testthat)

# Helper: build a minimal peptable with one modified and one unmodified row
# for the same sequence, plus quantification columns
make_peptable <- function(mod_vals, unmod_vals, ptmpos = list(c(3L)),
                          ptmtype = list(c("ph")),
                          quant_cols = c("S1", "S2", "S3")) {
  n_cols <- length(quant_cols)
  stopifnot(length(mod_vals) == n_cols, length(unmod_vals) == n_cols)

  mod_row <- data.frame(
    Sequence = "AASPEPR",
    stringsAsFactors = FALSE
  )
  mod_row[quant_cols] <- as.list(mod_vals)
  mod_row$PTMType <- ptmtype
  mod_row$PTMPos  <- ptmpos
  mod_row$Accession <- list("P12345")

  unmod_row <- data.frame(
    Sequence = "AASPEPR",
    stringsAsFactors = FALSE
  )
  unmod_row[quant_cols] <- as.list(unmod_vals)
  unmod_row$PTMType <- list(character(0))
  unmod_row$PTMPos  <- list(integer(0))
  unmod_row$Accession <- list("P12345")

  rbind(mod_row, unmod_row)
}

# Parameters helper
make_params <- function(quant_cols = c("S1", "S2", "S3")) {
  list(QuantColnames = quant_cols)
}

# ──────────────────────────────────────────────────────────────────────────────
# Basic correctness
# ──────────────────────────────────────────────────────────────────────────────

test_that("calcPTMOccupancy returns one row per modified peptidoform", {
  pep <- make_peptable(mod_vals = c(1, 2, 3), unmod_vals = c(1, 2, 3))
  occ <- calcPTMOccupancy(pep, make_params())

  expect_equal(nrow(occ), 1L)
})

test_that("occupancy is 0.5 when modified and unmodified intensities are equal", {
  # log2(x) == log2(x)  =>  x/(x+x) = 0.5
  pep <- make_peptable(mod_vals = c(2, 2, 2), unmod_vals = c(2, 2, 2))
  occ <- calcPTMOccupancy(pep, make_params())

  expect_equal(as.numeric(occ[1, c("S1", "S2", "S3")]),
               c(0.5, 0.5, 0.5), tolerance = 1e-9)
})

test_that("occupancy approaches 1 when modified >> unmodified", {
  # mod = 2^10 = 1024, unmod = 2^0 = 1  =>  1024/1025 ≈ 0.999
  pep <- make_peptable(mod_vals = c(10, 10, 10), unmod_vals = c(0, 0, 0))
  occ <- calcPTMOccupancy(pep, make_params())

  expect_true(all(as.numeric(occ[1, c("S1", "S2", "S3")]) > 0.99))
})

test_that("occupancy approaches 0 when modified << unmodified", {
  # mod = 2^0 = 1, unmod = 2^10 = 1024  =>  1/1025 ≈ 0.001
  pep <- make_peptable(mod_vals = c(0, 0, 0), unmod_vals = c(10, 10, 10))
  occ <- calcPTMOccupancy(pep, make_params())

  expect_true(all(as.numeric(occ[1, c("S1", "S2", "S3")]) < 0.01))
})

# ──────────────────────────────────────────────────────────────────────────────
# Output structure
# ──────────────────────────────────────────────────────────────────────────────

test_that("output contains expected columns", {
  pep <- make_peptable(mod_vals = c(1, 2, 3), unmod_vals = c(1, 2, 3))
  occ <- calcPTMOccupancy(pep, make_params())

  expect_true(all(c("Sequence", "Accession", "PTMPos", "PTMType", "S1", "S2", "S3") %in% names(occ)))
})

test_that("occupancy values are in [0, 1]", {
  pep <- make_peptable(mod_vals = c(1, 3, 5), unmod_vals = c(2, 4, 6))
  occ <- calcPTMOccupancy(pep, make_params())

  vals <- as.numeric(occ[1, c("S1", "S2", "S3")])
  expect_true(all(!is.na(vals)))
  expect_true(all(vals >= 0 & vals <= 1))
})

# ──────────────────────────────────────────────────────────────────────────────
# Edge cases
# ──────────────────────────────────────────────────────────────────────────────

test_that("returns empty data.frame when no modified peptides exist", {
  unmod_row <- data.frame(Sequence = "PEPTIDE", stringsAsFactors = FALSE)
  unmod_row[c("S1", "S2")] <- list(1, 2)
  unmod_row$PTMType <- list(character(0))
  unmod_row$PTMPos  <- list(integer(0))
  unmod_row$Accession <- list("P99999")

  occ <- calcPTMOccupancy(unmod_row, list(QuantColnames = c("S1", "S2")))
  expect_equal(nrow(occ), 0L)
})

test_that("returns empty data.frame when modified peptide has no unmodified counterpart", {
  mod_row <- data.frame(Sequence = "UNIQUEPEP", stringsAsFactors = FALSE)
  mod_row[c("S1", "S2")] <- list(3, 4)
  mod_row$PTMType <- list("ph")
  mod_row$PTMPos  <- list(2L)
  mod_row$Accession <- list("P00000")

  occ <- calcPTMOccupancy(mod_row, list(QuantColnames = c("S1", "S2")))
  expect_equal(nrow(occ), 0L)
})

test_that("NA in modified sample propagates to occupancy", {
  pep <- make_peptable(mod_vals = c(NA, 2, 3), unmod_vals = c(1, 2, 3))
  occ <- calcPTMOccupancy(pep, make_params())

  expect_true(is.na(occ[1, "S1"]))
  expect_false(is.na(occ[1, "S2"]))
})

test_that("NA in unmodified sample propagates to occupancy", {
  pep <- make_peptable(mod_vals = c(1, 2, 3), unmod_vals = c(NA, 2, 3))
  occ <- calcPTMOccupancy(pep, make_params())

  expect_true(is.na(occ[1, "S1"]))
  expect_false(is.na(occ[1, "S2"]))
})

# ──────────────────────────────────────────────────────────────────────────────
# Multiple unmodified rows: their signals are averaged
# ──────────────────────────────────────────────────────────────────────────────

test_that("multiple unmodified rows are averaged before occupancy calculation", {
  quant_cols <- c("S1")
  # Two unmodified rows with log2 intensities 0 and 2:
  # linear: 2^0=1 and 2^2=4  =>  column mean = (1+4)/2 = 2.5
  # modified row is set to log2(2.5), giving linear intensity exactly 2.5
  # occupancy = 2.5 / (2.5 + 2.5) = 0.5
  mod_row <- data.frame(Sequence = "MYPEPTIDE", stringsAsFactors = FALSE)
  mod_row[quant_cols] <- log2(2.5)
  mod_row$PTMType <- list("ph")
  mod_row$PTMPos  <- list(1L)
  mod_row$Accession <- list("P11111")

  unmod1 <- data.frame(Sequence = "MYPEPTIDE", stringsAsFactors = FALSE)
  unmod1[quant_cols] <- 0   # 2^0 = 1
  unmod1$PTMType <- list(character(0))
  unmod1$PTMPos  <- list(integer(0))
  unmod1$Accession <- list("P11111")

  unmod2 <- data.frame(Sequence = "MYPEPTIDE", stringsAsFactors = FALSE)
  unmod2[quant_cols] <- 2   # 2^2 = 4
  unmod2$PTMType <- list(character(0))
  unmod2$PTMPos  <- list(integer(0))
  unmod2$Accession <- list("P11111")

  pep <- rbind(mod_row, unmod1, unmod2)
  occ <- calcPTMOccupancy(pep, list(QuantColnames = quant_cols))

  expect_equal(nrow(occ), 1L)
  expect_equal(as.numeric(occ[1, "S1"]), 0.5, tolerance = 1e-9)
})

test_that("NA in any unmodified row propagates to occupancy for that sample", {
  quant_cols <- c("S1", "S2")

  mod_row <- data.frame(Sequence = "PEPTIDE", stringsAsFactors = FALSE)
  mod_row[quant_cols] <- list(2, 2)
  mod_row$PTMType <- list("ph")
  mod_row$PTMPos  <- list(1L)
  mod_row$Accession <- list("P33333")

  # Two unmodified rows: second has NA in S1
  unmod1 <- data.frame(Sequence = "PEPTIDE", stringsAsFactors = FALSE)
  unmod1[quant_cols] <- list(2, 2)
  unmod1$PTMType <- list(character(0))
  unmod1$PTMPos  <- list(integer(0))
  unmod1$Accession <- list("P33333")

  unmod2 <- data.frame(Sequence = "PEPTIDE", stringsAsFactors = FALSE)
  unmod2[quant_cols] <- list(NA_real_, 2)
  unmod2$PTMType <- list(character(0))
  unmod2$PTMPos  <- list(integer(0))
  unmod2$Accession <- list("P33333")

  pep <- rbind(mod_row, unmod1, unmod2)
  occ <- calcPTMOccupancy(pep, list(QuantColnames = quant_cols))

  # S1 has NA in one unmodified row -> occupancy must be NA for S1
  expect_true(is.na(occ[1, "S1"]))
  # S2 has no NA -> occupancy is valid
  expect_false(is.na(occ[1, "S2"]))
})

# ──────────────────────────────────────────────────────────────────────────────
# Multiple modified peptidoforms for the same sequence
# ──────────────────────────────────────────────────────────────────────────────

test_that("multiple modified rows for the same sequence each get their own row", {
  quant_cols <- c("S1")

  mod1 <- data.frame(Sequence = "STYPEP", stringsAsFactors = FALSE)
  mod1[quant_cols] <- 2
  mod1$PTMType <- list("ph")
  mod1$PTMPos  <- list(1L)
  mod1$Accession <- list("P22222")

  mod2 <- data.frame(Sequence = "STYPEP", stringsAsFactors = FALSE)
  mod2[quant_cols] <- 3
  mod2$PTMType <- list("ph")
  mod2$PTMPos  <- list(2L)
  mod2$Accession <- list("P22222")

  unmod <- data.frame(Sequence = "STYPEP", stringsAsFactors = FALSE)
  unmod[quant_cols] <- 2
  unmod$PTMType <- list(character(0))
  unmod$PTMPos  <- list(integer(0))
  unmod$Accession <- list("P22222")

  pep <- rbind(mod1, mod2, unmod)
  occ <- calcPTMOccupancy(pep, list(QuantColnames = quant_cols))

  expect_equal(nrow(occ), 2L)
})
