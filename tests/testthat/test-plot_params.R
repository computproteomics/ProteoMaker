library(testthat)

# Mock the param_table function for testing
param_table <- function() {
    data.frame(Group = c("paramGroundTruth", "paramGroundTruth", "paramProteoformAb", "paramDigest"),
               row.names = c("NumReps", "NumCond", "ProteoformAb", "DigestionEnzyme"))
}

test_that("plot_params returns a plotly object", {
    BenchMatrix <- data.frame(NumReps = c(2, 4), NumCond = c(3, 2))
    result <- plot_params(BenchMatrix, current_row = 1)
    
    expect_true(inherits(result, "plotly"))
})


test_that("plot_params handles rows with NA values", {
    BenchMatrix <- data.frame(NumReps = c(2, NA), NumCond = c(3, NA))
    result <- plot_params(BenchMatrix, current_row = 1)
    
    expect_true(inherits(result, "plotly"))
})
