library(testthat)
test_that("generate_combinations returns correct combinations for single parameter", {
    params <- list(NumReps = c(2, 4, 6), NumCond = c(2,3))
    combinations <- generate_combinations(params)
    
    expect_equal(length(combinations), 6)
    expect_equal(combinations[[1]]$NumReps, 2)
    expect_equal(combinations[[2]]$NumReps, 4)
    expect_equal(combinations[[3]]$NumReps, 6)
})

test_that("generate_combinations returns correct combinations for multiple parameters", {
    params <- list(NumReps = c(2, 4), NumCond = c(1, 2))
    combinations <- generate_combinations(params)
    
    expect_equal(length(combinations), 4)
    expect_equal(combinations[[1]]$NumReps, 2)
    expect_equal(combinations[[1]]$NumCond, 1)
    expect_equal(combinations[[2]]$NumReps, 4)
    expect_equal(combinations[[2]]$NumCond, 1)
    expect_equal(combinations[[3]]$NumReps, 2)
    expect_equal(combinations[[3]]$NumCond, 2)
    expect_equal(combinations[[4]]$NumReps, 4)
    expect_equal(combinations[[4]]$NumCond, 2)
})

