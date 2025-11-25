library(testthat)


test_that("gather_all_sims returns correct structure", {
  Param <- def_param()
  Param$paramGroundTruth$NumReps <- 5
  ll <- list.files(tempdir(), pattern = "output", full.names = TRUE)
  unlink(ll, recursive = TRUE)
  proteomaker_config <- set_proteomaker(resultFilePath = tempdir())
  benchmarks <- run_sims(Param, proteomaker_config)

  sims <- gather_all_sims(proteomaker_config)

  expect_true(length(sims[[1]]$Benchmarks) > 0)
  expect_true(length(sims[[1]]$Param) > 0)

  sims <- gather_all_sims(proteomaker_config, stage = "MSRun")
  expect_true(length(sims[[1]]$Benchmarks) == 0)
  expect_true(length(sims[[1]]$Param) > 0)
})

test_that("gather_all_sims handles empty directories", {
  ll <- list.files(tempdir(), pattern = "output", full.names = TRUE)
  unlink(ll, recursive = TRUE)
  # Mock configuration with a different temp directory (empty)
  proteomaker_config <- set_proteomaker(resultFilePath = tempdir())

  # Run the function on an empty directory
  all_results_empty <- gather_all_sims(proteomaker_config, stage = "DataAnalysis")

  # Check that the result is an empty list
  expect_equal(length(all_results_empty), 0)
})
