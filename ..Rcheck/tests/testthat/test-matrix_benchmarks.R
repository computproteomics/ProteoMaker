library(testthat)

test_that("matrix_benchmarks returns a correctly formatted data frame", {
  Param <- def_param()
  Param$paramGroundTruth$NumReps <- 5
  ll <- list.files(tempdir(), pattern = "output", full.names = TRUE)
  unlink(ll, recursive = TRUE)
  proteomaker_config <- test_proteomaker_config(resultFilePath = tempdir())
  benchmarks <- run_sims(Param, proteomaker_config)
  results <- matrix_benchmarks(benchmarks, proteomaker_config)

  expect_true(is.data.frame(results))
  expect_true(results$NumReps == 5)
})
