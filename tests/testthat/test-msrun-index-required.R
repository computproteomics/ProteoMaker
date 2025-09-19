test_that("MSRunSim requires non-NULL searchIndex", {
  # Minimal digested table with required columns
  dig <- data.frame(
    Sequence = c("PEPTIDEA", "PEPTIDEB"),
    Accession = I(list("P1", "P2")),
    PTMPos = I(list(integer(0), integer(0))),
    PTMType = I(list(character(0), character(0))),
    C_1_R_1 = c(10, 11),
    C_2_R_1 = c(12, 13),
    stringsAsFactors = FALSE
  )
  params <- list(
    MSNoise = 0,
    PercDetectability = 1,
    PercDetectedVal = 1,
    WeightDetectVal = 0,
    WrongIDs = 0.1,
    QuantColnames = c("C_1_R_1", "C_2_R_1")
  )
  expect_error(
    MSRunSim(dig, params, searchIndex = NULL),
    regexp = "requires a non-NULL searchIndex"
  )
})

