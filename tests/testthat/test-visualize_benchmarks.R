library(testthat)

test_that("visualize_benchmarks returns a plotly object", {
  Param <- def_param()
  ll <- list.files(tempdir(), pattern = "output", full.names = TRUE)
  unlink(ll, recursive = TRUE)
  proteomaker_config <- test_proteomaker_config(resultFilePath = tempdir())
  benchmarks <- run_sims(Param, proteomaker_config)
  benchmatrix <- matrix_benchmarks(benchmarks, proteomaker_config)
  result <- visualize_benchmarks(benchmatrix)
})

test_that("visualize_benchmarks supports shared compare_par colors", {
  benchmatrix <- data.frame(
    WrongIDs = c(0.01, 0.01, 0.01, 0.01, 0.05, 0.05, 0.05, 0.05),
    NumReps = c(2, 2, 3, 3, 2, 2, 3, 3),
    numPeptides = c(100, 105, 120, 125, 90, 95, 110, 115)
  )

  expect_no_error(
    visualize_benchmarks(
      benchmatrix,
      benchmarks = "numPeptides",
      ref_par = "WrongIDs",
      compare_par = "NumReps",
      errorbar = TRUE
    )
  )

  expect_no_error(
    visualize_benchmarks(
      benchmatrix,
      benchmarks = "numPeptides",
      ref_par = "WrongIDs",
      compare_par = "NumReps",
      errorbar = TRUE,
      errorstyle = "area"
    )
  )
})

test_that("visualize_benchmarks supports area uncertainty without compare_par", {
  benchmatrix <- data.frame(
    WrongIDs = c(0.01, 0.01, 0.05, 0.05),
    numPeptides = c(100, 105, 90, 95)
  )

  expect_no_error(
    visualize_benchmarks(
      benchmatrix,
      benchmarks = "numPeptides",
      ref_par = "WrongIDs",
      errorbar = TRUE,
      errorstyle = "area"
    )
  )
})
