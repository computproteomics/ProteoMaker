test_that("calcBenchmarks handles FDR columns and small adjusted sets", {
  Param <- list(
    NumCond = 2,
    NumReps = 1,
    QuantColnames = c("C_1_R_1", "C_2_R_1")
  )
  # Minimal protein stats with one FDR column and a log-ratio
  Stats <- data.frame(
    `log-ratios 2 vs 1` = c(0.5, -0.2, 0.0),
    `FDR_limma Condition 2 vs Condition 1` = c(0.001, 0.5, NA_real_),
    min1Reg = c(TRUE, FALSE, FALSE),
    Regulation_Pattern = c("1;1", "1;1", "1;1"),
    Regulation_Amplitude = c("NA;NA", "NA;NA", "NA;NA"),
    Sequence = c("PEPA;PEPB", "PEPC;PEPD", "PEPE;PEPF"),
    num_accs = c("1", "1", "1"),
    WrongID = c("FALSE", "FALSE", "FALSE"),
    Proteoform_ID = c("1", "1", "1"),
    MC = c("0", "0", "0"),
    C_1_R_1 = c(10, 11, 12),
    C_2_R_1 = c(10.5, 10.9, 12.1),
    stringsAsFactors = FALSE
  )
  # Minimal peptide stats with matching columns
  StatsPep <- data.frame(
    Sequence = c("PEPA", "PEPB", "PEPC"),
    Accession = I(list("P1", "P2", "P3")),
    MC = I(list(0L, 0L, 0L)),
    PTMPos = I(list(integer(0), integer(0), integer(0))),
    PTMType = I(list(character(0), character(0), character(0))),
    min1Reg = c(TRUE, FALSE, FALSE),
    `log-ratios 2 vs 1` = c(0.4, -0.1, 0.0),
    `FDR_limma Condition 2 vs Condition 1` = c(0.002, 0.6, NA_real_),
    Regulation_Pattern = c("1;1", "1;1", "1;1"),
    Regulation_Amplitude = c("0.5;0.5", "0;0", "NA;NA"),
    C_1_R_1 = c(8, 7, 9),
    C_2_R_1 = c(8.4, 7.1, 9.0),
    stringsAsFactors = FALSE
  )
  bm <- calcBenchmarks(Stats, StatsPep, Param)
  expect_type(bm, "list")
  expect_true("globalBMs" %in% names(bm))
})

test_that("calcBenchmarks handles NumCond > 2 without vector-replacement warnings", {
  Param <- list(
    NumCond = 3,
    NumReps = 1,
    QuantColnames = c("C_1_R_1", "C_2_R_1", "C_3_R_1")
  )
  Stats <- data.frame(
    `log-ratios 2 vs 1` = c(1.0, -1.0, 0.0),
    `log-ratios 3 vs 1` = c(2.0, -2.0, 0.0),
    `FDR_limma Condition 2 vs Condition 1` = c(0.001, 0.5, NA_real_),
    `FDR_limma Condition 3 vs Condition 1` = c(0.001, 0.5, NA_real_),
    min1Reg = c(TRUE, FALSE, FALSE),
    allReg = c(TRUE, FALSE, FALSE),
    Regulation_Pattern = c("1;1;1", "1;1;1", "1;1;1"),
    Regulation_Amplitude = c("1;NA", "1;NA", "NA;NA"),
    Sequence = c("PEPA;PEPB", "PEPC;PEPD", "PEPE;PEPF"),
    num_accs = c("1", "1", "1"),
    WrongID = c("FALSE", "FALSE", "FALSE"),
    Proteoform_ID = c("1", "1", "1"),
    MC = c("0", "0", "0"),
    C_1_R_1 = c(10, 11, 12),
    C_2_R_1 = c(11.0, 10.0, 12.0),
    C_3_R_1 = c(12.0, 9.0, 12.0),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  StatsPep <- data.frame(
    Sequence = c("PEPA", "PEPB", "PEPC"),
    Accession = I(list("P1", "P2", "P3")),
    MC = I(list(0L, 0L, 0L)),
    PTMPos = I(list(integer(0), integer(0), integer(0))),
    PTMType = I(list(character(0), character(0), character(0))),
    min1Reg = c(TRUE, FALSE, FALSE),
    `log-ratios 2 vs 1` = c(1.0, -1.0, 0.0),
    `log-ratios 3 vs 1` = c(2.0, -2.0, 0.0),
    `FDR_limma Condition 2 vs Condition 1` = c(0.002, 0.6, NA_real_),
    `FDR_limma Condition 3 vs Condition 1` = c(0.002, 0.6, NA_real_),
    Regulation_Pattern = c("1;1;1", "1;1;1", "1;1;1"),
    Regulation_Amplitude = c("1;NA", "1;NA", "NA;NA"),
    C_1_R_1 = c(8, 7, 9),
    C_2_R_1 = c(9.0, 6.0, 9.0),
    C_3_R_1 = c(10.0, 5.0, 9.0),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  expect_no_warning(
    bm <- calcBenchmarks(Stats, StatsPep, Param)
  )
  expect_type(bm, "list")
  expect_true("globalBMs" %in% names(bm))
  expect_length(bm$globalBMs$meanSquareDiffFCProt, 1)
  expect_length(bm$globalBMs$meanSquareDiffFCPep, 1)
})
