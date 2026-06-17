library(testthat)

# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────

# Build a minimal proteoformAb table.
# rows: list of lists, each with: seq, acc, ptmpos, ptmtype, vals (named by quant cols)
make_proteoform_table <- function(rows, quant_cols) {
  frames <- lapply(rows, function(r) {
    df <- data.frame(Sequence = r$seq, Accession = r$acc,
                     stringsAsFactors = FALSE)
    df[quant_cols] <- as.list(r$vals)
    df$PTMPos  <- I(list(r$ptmpos))
    df$PTMType <- I(list(r$ptmtype))
    df
  })
  do.call(rbind, frames)
}

# ──────────────────────────────────────────────────────────────────────────────
# Basic correctness
# ──────────────────────────────────────────────────────────────────────────────

test_that("calcGroundTruthPTMOccupancy calculates replicate-level occupancy", {
  proteoforms <- data.frame(
    Sequence = rep("MSTK", 3),
    Accession = rep("P1", 3),
    C_1_R_1 = log2(c(25, 75, 50)),
    C_1_R_2 = log2(c(50, 50, 100)),
    stringsAsFactors = FALSE
  )
  proteoforms$PTMPos <- I(list(integer(0), 2L, 2L))
  proteoforms$PTMType <- I(list(character(0), "Phospho", "Phospho"))

  occ <- calcGroundTruthPTMOccupancy(
    proteoforms,
    list(QuantColnames = c("C_1_R_1", "C_1_R_2"))
  )

  expect_equal(nrow(occ), 1)
  expect_equal(occ$C_1_R_1, 125 / 150)
  expect_equal(occ$C_1_R_2, 150 / 200)
  expect_identical(occ$Accession[[1]], "P1")
  expect_identical(occ$PTMPos[[1]], 2L)
  expect_identical(occ$PTMType[[1]], "Phospho")
})

test_that("calcGroundTruthPTMOccupancy returns one row per PTM site", {
  proteoforms <- data.frame(
    Sequence = rep("MSTKST", 3),
    Accession = rep("P1", 3),
    C_1_R_1 = log2(c(25, 25, 50)),
    stringsAsFactors = FALSE
  )
  proteoforms$PTMPos <- I(list(integer(0), c(2L, 5L), c(5L, 2L)))
  proteoforms$PTMType <- I(list(character(0), c("Phospho", "Phospho"), c("Phospho", "Phospho")))

  occ <- calcGroundTruthPTMOccupancy(
    proteoforms,
    list(QuantColnames = "C_1_R_1")
  )

  expect_equal(nrow(occ), 2)
  expect_equal(occ$C_1_R_1, c(75 / 100, 75 / 100))
  expect_equal(sort(unlist(occ$PTMPos)), c(2L, 5L))
  expect_true(all(unlist(occ$PTMType) == "Phospho"))
})

# ──────────────────────────────────────────────────────────────────────────────
# Empty / no-PTM inputs
# ──────────────────────────────────────────────────────────────────────────────

test_that("returns empty data.frame for NULL input", {
  occ <- calcGroundTruthPTMOccupancy(NULL, list(QuantColnames = "C_1_R_1"))
  expect_equal(nrow(occ), 0L)
})

test_that("returns empty data.frame for zero-row input", {
  empty <- data.frame(Sequence = character(0), Accession = character(0),
                      C_1_R_1 = numeric(0), stringsAsFactors = FALSE)
  empty$PTMPos  <- I(list())
  empty$PTMType <- I(list())
  occ <- calcGroundTruthPTMOccupancy(empty, list(QuantColnames = "C_1_R_1"))
  expect_equal(nrow(occ), 0L)
})

test_that("returns empty data.frame when no proteoforms are modified", {
  # All PTMType entries are empty — function should message and return empty
  pf <- make_proteoform_table(
    list(
      list(seq = "AAPK", acc = "P1", ptmpos = integer(0), ptmtype = character(0), vals = c(C_1_R_1 = 4)),
      list(seq = "SEQB", acc = "P1", ptmpos = integer(0), ptmtype = character(0), vals = c(C_1_R_1 = 4))
    ),
    "C_1_R_1"
  )
  occ <- calcGroundTruthPTMOccupancy(pf, list(QuantColnames = "C_1_R_1"))
  expect_equal(nrow(occ), 0L)
})

# ──────────────────────────────────────────────────────────────────────────────
# Missing values
# ──────────────────────────────────────────────────────────────────────────────

test_that("returns NA occupancy when all denominator values are NA", {
  # All three proteoforms for P1 have NA in C_1_R_1 → denominator is all-NA
  pf <- make_proteoform_table(
    list(
      list(seq = "MSTK", acc = "P1", ptmpos = integer(0),  ptmtype = character(0), vals = c(C_1_R_1 = NA_real_)),
      list(seq = "MSTK", acc = "P1", ptmpos = 2L,          ptmtype = "Phospho",    vals = c(C_1_R_1 = NA_real_))
    ),
    "C_1_R_1"
  )
  occ <- calcGroundTruthPTMOccupancy(pf, list(QuantColnames = "C_1_R_1"))
  expect_equal(nrow(occ), 1L)
  expect_true(is.na(occ$C_1_R_1))
})

test_that("NA in one replicate does not affect the other replicate's occupancy", {
  # C_1_R_1: unmod=NA, mod=log2(75) → denominator only has mod row → occ = 75/75 = 1
  # C_1_R_2: unmod=log2(25), mod=log2(75) → occ = 75/100 = 0.75
  pf <- make_proteoform_table(
    list(
      list(seq = "MSTK", acc = "P1", ptmpos = integer(0), ptmtype = character(0),
           vals = c(C_1_R_1 = NA_real_, C_1_R_2 = log2(25))),
      list(seq = "MSTK", acc = "P1", ptmpos = 2L,         ptmtype = "Phospho",
           vals = c(C_1_R_1 = log2(75), C_1_R_2 = log2(75)))
    ),
    c("C_1_R_1", "C_1_R_2")
  )
  occ <- calcGroundTruthPTMOccupancy(pf, list(QuantColnames = c("C_1_R_1", "C_1_R_2")))
  expect_equal(nrow(occ), 1L)
  expect_equal(occ$C_1_R_1, 1,    tolerance = 1e-9)
  expect_equal(occ$C_1_R_2, 0.75, tolerance = 1e-9)
})

# ──────────────────────────────────────────────────────────────────────────────
# Multiple proteins — accessions must not contaminate each other
# ──────────────────────────────────────────────────────────────────────────────

test_that("occupancy is computed independently per protein", {
  # Protein P_AAA: unmod=log2(50), mod=log2(50) → occ = 50/100 = 0.5
  # Protein P_BBB: unmod=log2(10), mod=log2(90) → occ = 90/100 = 0.9
  # If accessions were mixed, the denominator would be 200 and numerators wrong.
  pf <- make_proteoform_table(
    list(
      list(seq = "SEQA", acc = "P_AAA", ptmpos = integer(0), ptmtype = character(0), vals = c(C_1_R_1 = log2(50))),
      list(seq = "SEQA", acc = "P_AAA", ptmpos = 2L,         ptmtype = "Phospho",    vals = c(C_1_R_1 = log2(50))),
      list(seq = "SEQB", acc = "P_BBB", ptmpos = integer(0), ptmtype = character(0), vals = c(C_1_R_1 = log2(10))),
      list(seq = "SEQB", acc = "P_BBB", ptmpos = 3L,         ptmtype = "Phospho",    vals = c(C_1_R_1 = log2(90)))
    ),
    "C_1_R_1"
  )
  occ <- calcGroundTruthPTMOccupancy(pf, list(QuantColnames = "C_1_R_1"))

  expect_equal(nrow(occ), 2L)
  occ_a <- occ[unlist(occ$Accession) == "P_AAA", ]
  occ_b <- occ[unlist(occ$Accession) == "P_BBB", ]
  expect_equal(occ_a$C_1_R_1, 0.5, tolerance = 1e-9)
  expect_equal(occ_b$C_1_R_1, 0.9, tolerance = 1e-9)
})

test_that("protein with no modified proteoforms contributes no rows", {
  # P_AAA has a modified proteoform; P_BBB has only unmodified ones.
  # Output should have exactly one row for P_AAA only.
  pf <- make_proteoform_table(
    list(
      list(seq = "SEQA", acc = "P_AAA", ptmpos = integer(0), ptmtype = character(0), vals = c(C_1_R_1 = log2(50))),
      list(seq = "SEQA", acc = "P_AAA", ptmpos = 2L,         ptmtype = "Phospho",    vals = c(C_1_R_1 = log2(50))),
      list(seq = "SEQB", acc = "P_BBB", ptmpos = integer(0), ptmtype = character(0), vals = c(C_1_R_1 = log2(80)))
    ),
    "C_1_R_1"
  )
  occ <- calcGroundTruthPTMOccupancy(pf, list(QuantColnames = "C_1_R_1"))

  expect_equal(nrow(occ), 1L)
  expect_identical(unlist(occ$Accession), "P_AAA")
})
