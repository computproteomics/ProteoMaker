library(testthat)


test_that("gather_all_sims returns correct structure", {
  Param <- def_param()
  Param$paramGroundTruth$NumReps <- 5
  ll <- list.files(tempdir(), pattern="output", full.names = TRUE)
  unlink(ll, recursive = TRUE)
  phosfake_config <- set_phosfake(resultFilePath = tempdir())
  benchmarks <- run_sims(Param, phosfake_config)

  sims <- gather_all_sims(phosfake_config)

  expect_true(length(sims[[1]]$Benchmarks) > 0)
  expect_true(length(sims[[1]]$Param) > 0)

  sims <- gather_all_sims(phosfake_config, stage = "MSRun")
  expect_true(length(sims[[1]]$Benchmarks) == 0)
  expect_true(length(sims[[1]]$Param) > 0)


})

test_that("gather_all_sims handles empty directories", {
  ll <- list.files(tempdir(), pattern="output", full.names = TRUE)
  unlink(ll, recursive = TRUE)
  # Mock configuration with a different temp directory (empty)
  phosfake_config <- set_phosfake(resultFilePath = tempdir())

  # Run the function on an empty directory
  all_results_empty <- gather_all_sims(phosfake_config, stage = "DataAnalysis")

  # Check that the result is an empty list
  expect_equal(length(all_results_empty), 0)
})


