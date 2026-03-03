test_that("matrix_benchmarks exposes base adjusted PTM metric names", {
  # minimal synthetic run-level object with globalBMs only
  allBs <- list(
    sim1 = list(
      Benchmarks = list(
        globalBMs = list(
          aucDiffRegAdjModPep = 0.81,
          tprAdjModPep0.01 = 0.42,
          tprAdjModPep0.05 = 0.57,
          tFDRAdjModPep0.01 = 0.03,
          tFDRAdjModPep0.05 = 0.06
        )
      ),
      Param = list(
        WrongIDs = 0.01,
        NumReps = 3
      )
    )
  )

  bm <- matrix_benchmarks(allBs, Config = list())

  expect_true(all(c(
    "aucDiffRegAdjModPep",
    "tprAdjModPep0.01",
    "tprAdjModPep0.05",
    "tFDRAdjModPep0.01",
    "tFDRAdjModPep0.05"
  ) %in% colnames(bm)))

  expect_equal(unname(bm[1, "aucDiffRegAdjModPep"]), 0.81)
  expect_equal(unname(bm[1, "tprAdjModPep0.01"]), 0.42)
  expect_equal(unname(bm[1, "tprAdjModPep0.05"]), 0.57)
  expect_equal(unname(bm[1, "tFDRAdjModPep0.01"]), 0.03)
  expect_equal(unname(bm[1, "tFDRAdjModPep0.05"]), 0.06)
})
