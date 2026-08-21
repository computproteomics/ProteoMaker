library(testthat)

test_that("visualize_benchmarks returns a plotly object", {
  benchmatrix <- data.frame(
    WrongIDs = c(0.01, 0.05),
    numPeptides = c(100, 90)
  )

  expect_no_error(
    visualize_benchmarks(
      benchmatrix,
      benchmarks = "numPeptides",
      ref_par = "WrongIDs"
    )
  )
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
