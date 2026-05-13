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
