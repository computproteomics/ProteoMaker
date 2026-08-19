library(testthat)

# Mock the param_table function for testing
param_table <- function() {
  data.frame(
    Category = c("Experimental Design", "Experimental Design", "Proteoform Quantities", "Enzymatic Digestion"),
    Group = c("paramGroundTruth", "paramGroundTruth", "paramProteoformAb", "paramDigest"),
    Explanation = c("Number of replicates.", "Number of conditions.", "Abundance setting.", "Digestion enzyme."),
    MinValue = c(1, 1, NA, NA),
    MaxValue = c(10, 10, NA, NA),
    DefaultValue = c(3, 2, NA, "trypsin"),
    row.names = c("NumReps", "NumCond", "ProteoformAb", "Enzyme")
  )
}

test_that("plot_params returns a plotly object", {
  BenchMatrix <- data.frame(NumReps = c(2, 4), NumCond = c(3, 2))
  result <- plot_params(BenchMatrix, current_row = 1)

  expect_true(inherits(result, "plotly"))
})

test_that("plot_params includes nonnumeric parameters as a table", {
  BenchMatrix <- data.frame(NumReps = 2, Enzyme = "trypsin")
  result <- plot_params(BenchMatrix, current_row = 1)

  expect_true(inherits(result, "plotly"))
  built <- plotly::plotly_build(result)
  expect_true(any(vapply(built$x$data, function(x) identical(x$type, "table"), logical(1))))
})

test_that("render_parameter_table returns selected parameter columns", {
  result <- render_parameter_table(
    list(paramGroundTruth = list(NumReps = 2), paramDigest = list(Enzyme = "trypsin")),
    print_table = FALSE
  )

  expect_named(result, c("Parameter", "Value", "Description"))
  expect_equal(result$Parameter, c("NumReps", "Enzyme"))
  expect_equal(result$Value, c("2", "trypsin"))
})


test_that("plot_params handles rows with NA values", {
  BenchMatrix <- data.frame(NumReps = c(2, NA), NumCond = c(3, NA))
  result <- plot_params(BenchMatrix, current_row = 1)

  expect_true(inherits(result, "plotly"))
})
