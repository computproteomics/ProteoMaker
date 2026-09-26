# Clean-data fold-change benchmarks: with zero noise the measured log-ratios must
# exactly match the ground-truth fold changes that calcBenchmarks derives from
# Regulation_Pattern/Regulation_Amplitude. Replicates are identical (zero variance),
# so one "decoy" row with within-condition replicate spread is included to give
# limma's eBayes finite residual variance; the decoy's own fold change is exact,
# hence it contributes zero benchmark error.
# Convention being pinned down: protein-level pattern cells drop unregulated
# proteoforms (upstream paste(unlist())), peptide-level cells keep them as zero rows.

Param2 <- list(
  NumCond = 2, NumReps = 2, StatPaired = FALSE,
  QuantColnames = c("C_1_R_1", "C_1_R_2", "C_2_R_1", "C_2_R_2")
)
Param3 <- list(
  NumCond = 3, NumReps = 2, StatPaired = FALSE,
  QuantColnames = c("C_1_R_1", "C_1_R_2", "C_2_R_1", "C_2_R_2", "C_3_R_1", "C_3_R_2")
)

fc_pep_row <- function(Sequence, c1, c2, pattern, amp, c3 = NULL, acc = "P1") {
  # c1/c2/c3: two replicate log2 intensities per condition
  row <- data.frame(
    Sequence = Sequence,
    Accession = I(list(acc)),
    MC = I(list(0L)),
    PTMPos = I(list(integer(0))),
    PTMType = I(list(character(0))),
    `log-ratios 2 vs 1` = NA_real_,
    `FDR_limma Condition 2 vs Condition 1` = NA_real_,
    Regulation_Pattern = pattern,
    Regulation_Amplitude = amp,
    C_1_R_1 = c1[1], C_1_R_2 = c1[2], C_2_R_1 = c2[1], C_2_R_2 = c2[2],
    stringsAsFactors = FALSE, check.names = FALSE
  )
  if (!is.null(c3)) {
    row <- cbind(row, data.frame(
      `log-ratios 3 vs 1` = NA_real_,
      `FDR_limma Condition 3 vs Condition 1` = NA_real_,
      C_3_R_1 = c3[1], C_3_R_2 = c3[2],
      stringsAsFactors = FALSE, check.names = FALSE
    ))
  }
  row
}

fc_prot_row <- function(c1, c2, pattern, amp, c3 = NULL) {
  row <- data.frame(
    `log-ratios 2 vs 1` = NA_real_,
    `FDR_limma Condition 2 vs Condition 1` = NA_real_,
    Regulation_Pattern = pattern,
    Regulation_Amplitude = amp,
    Sequence = "PEPA",
    num_accs = "1",
    WrongID = "FALSE",
    Proteoform_ID = "1",
    MC = "0",
    C_1_R_1 = c1[1], C_1_R_2 = c1[2], C_2_R_1 = c2[1], C_2_R_2 = c2[2],
    stringsAsFactors = FALSE, check.names = FALSE
  )
  if (!is.null(c3)) {
    row <- cbind(row, data.frame(
      `log-ratios 3 vs 1` = NA_real_,
      `FDR_limma Condition 3 vs Condition 1` = NA_real_,
      C_3_R_1 = c3[1], C_3_R_2 = c3[2],
      stringsAsFactors = FALSE, check.names = FALSE
    ))
  }
  row
}

# decoys: shifted by exactly 1 log2 unit per comparison, with replicate spread
decoy_prot2 <- fc_prot_row(c(10, 10.5), c(11, 11.5), "-0.5;0.5", "1")
decoy_prot3 <- fc_prot_row(c(10, 10.5), c(11, 11.5), "-1;0;1", "1", c3 = c(12, 12.5))
decoy_pep2 <- fc_pep_row("SDECOY", c(10, 10.5), c(11, 11.5), "c(-0.5, 0.5)", "1")

StatsPep_empty2 <- fc_pep_row("SDUMMY", c(10, 10), c(10, 10), "NULL", "NA")

test_that("zero-noise single-proteoform protein matches ground truth exactly", {
  # pattern -0.5;0.5, amplitude 1: condition 2 sits 1 log2 unit above condition 1
  Stats <- runPolySTest(rbind(decoy_prot2, fc_prot_row(c(8.5, 8.5), c(9.5, 9.5), "-0.5;0.5", "1")),
                        Param2, refCond = 1, onlyLIMMA = TRUE)
  bm <- calcBenchmarks(Stats, runPolySTest(StatsPep_empty2, Param2, refCond = 1, onlyLIMMA = TRUE), Param2)
  expect_equal(bm$globalBMs$meanSquareDiffFCProt, 0, tolerance = 1e-10)
})

test_that("zero-noise unregulated co-proteoform dilutes the group fold change exactly", {
  # regulated proteoform: 8 -> 10 (truth 2); unregulated co-proteoform stays at 9;
  # group mean: 8.5 -> 9.5 (measured 1); benchmark truth counts the regulated
  # proteoform only, so the dilution shows up as exactly (2 - 1)^2
  Stats <- runPolySTest(rbind(decoy_prot2, fc_prot_row(c(8.5, 8.5), c(9.5, 9.5), "-0.5;0.5", "2")),
                        Param2, refCond = 1, onlyLIMMA = TRUE)
  bm <- calcBenchmarks(Stats, runPolySTest(StatsPep_empty2, Param2, refCond = 1, onlyLIMMA = TRUE), Param2)
  expect_equal(bm$globalBMs$meanSquareDiffFCProt, ((2 - 1)^2 + 0) / 2, tolerance = 1e-10)
})

test_that("zero-noise unregulated protein group is excluded from the FC benchmark", {
  Stats <- runPolySTest(rbind(
    decoy_prot2,
    fc_prot_row(c(8.5, 8.5), c(9.5, 9.5), "-0.5;0.5", "1"),
    fc_prot_row(c(11, 11), c(11, 11), "", "NA")
  ), Param2, refCond = 1, onlyLIMMA = TRUE)
  bm <- calcBenchmarks(Stats, runPolySTest(StatsPep_empty2, Param2, refCond = 1, onlyLIMMA = TRUE), Param2)
  expect_equal(bm$globalBMs$meanSquareDiffFCProt, 0, tolerance = 1e-10)
})

test_that("zero-noise protein enters both comparisons at NumCond = 3", {
  # pattern -1;0;1: truth diffs are 1 (2 vs 1) and 2 (3 vs 1); measured are 0 and 1
  Stats <- runPolySTest(rbind(decoy_prot3, fc_prot_row(c(8, 8), c(8, 8), "-1;0;1", "1", c3 = c(9, 9))),
                        Param3, refCond = 1, onlyLIMMA = TRUE)
  bm <- calcBenchmarks(Stats, runPolySTest(rbind(
    fc_pep_row("SDECOY3", c(10, 10.5), c(11, 11.5), "c(-1, 0, 1)", "1", c3 = c(12, 12.5)),
    fc_pep_row("SDUMMY", c(10, 10), c(10, 10), "NULL", "NA", c3 = c(10, 10))
  ), Param3, refCond = 1, onlyLIMMA = TRUE), Param3)
  # clean row: truth (1, 2) vs measured (0, 1) -> mean of squared diffs = 1; decoy exact -> 0
  expect_equal(bm$globalBMs$meanSquareDiffFCProt, (1 + 0) / 2, tolerance = 1e-10)
})

test_that("zero-noise peptides match ground truth, including unregulated co-proteoforms", {
  # S1: single regulated proteoform, truth -0.4, measured -0.4 -> 0 error
  # S2: regulated proteoform plus "NULL" co-proteoform: truth mean(-0.6, 0) = -0.3,
  #     measured -0.3 -> 0 error
  # S3: fully unregulated, skipped
  StatsPep <- runPolySTest(rbind(
    decoy_pep2,
    fc_pep_row("S1", c(7.5, 7.5), c(7.1, 7.1), "c(0.5, -0.5)", "0.4", acc = "P1"),
    fc_pep_row("S2", c(8, 8), c(7.7, 7.7), "c(0.5, -0.5);NULL", "0.6", acc = "P2"),
    fc_pep_row("S3", c(10, 10), c(10, 10), "NULL", "NA", acc = "P3")
  ), Param2, refCond = 1, onlyLIMMA = TRUE)
  Stats <- runPolySTest(rbind(decoy_prot2, fc_prot_row(c(8.5, 8.5), c(9.5, 9.5), "-0.5;0.5", "1")),
                        Param2, refCond = 1, onlyLIMMA = TRUE)
  bm <- calcBenchmarks(Stats, StatsPep, Param2)
  expect_equal(bm$globalBMs$meanSquareDiffFCPep, 0, tolerance = 1e-10)

  # perturb S1's measurement by 0.2: benchmark must see exactly (-0.2 - (-0.4))^2 / 3
  # (S2 and the decoy still exact), pinning the truth average over the NULL co-proteoform
  StatsPep$C_2_R_1[StatsPep$Sequence == "S1"] <- 7.3
  StatsPep$C_2_R_2[StatsPep$Sequence == "S1"] <- 7.3
  bm2 <- calcBenchmarks(Stats, runPolySTest(StatsPep, Param2, refCond = 1, onlyLIMMA = TRUE), Param2)
  expect_equal(bm2$globalBMs$meanSquareDiffFCPep, 0.2^2 / 3, tolerance = 1e-10)
})
