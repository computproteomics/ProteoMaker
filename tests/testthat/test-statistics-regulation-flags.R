test_that("runPolySTest regulation flags distinguish any vs all regulated contributors", {
  params <- list(
    QuantColnames = c("C_1_R_1", "C_1_R_2", "C_2_R_1", "C_2_R_2"),
    NumCond = 2,
    NumReps = 2,
    StatPaired = FALSE
  )

  dat <- data.frame(
    C_1_R_1 = c(10, 10, 10, 10),
    C_1_R_2 = c(10, 10, 10, 10),
    C_2_R_1 = c(11, 11, 11, 11),
    C_2_R_2 = c(11, 11, 11, 11),
    stringsAsFactors = FALSE
  )
  dat$Regulation_Amplitude <- I(list(
    c(1, 2),
    c(1, NA),
    c(1, 0),
    c(NA, NA)
  ))
  dat$Regulation_Pattern <- I(rep(list(c(0.5, -0.5)), 4))

  stats <- runPolySTest(dat, params, refCond = 1, onlyLIMMA = TRUE, cores = 1)

  expect_identical(stats$min1Reg, c(TRUE, TRUE, TRUE, FALSE))
  expect_identical(stats$allReg, c(TRUE, FALSE, FALSE, FALSE))
})
