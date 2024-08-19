library(testthat)

test_that("visualize_benchmarks returns a plotly object", {
    Param <- def_param()
    ll <- list.files(tempdir(), pattern="output", full.names = TRUE)
    unlink(ll, recursive = TRUE)
    phosfake_config <- set_phosfake(resultFilePath = tempdir())
    benchmarks <- run_sims(Param, phosfake_config)
    benchmatrix <- matrix_benchmarks(benchmarks, phosfake_config)
    result <- suppressWarnings(visualize_benchmarks(benchmatrix, current_row = 1))
    expect_true(inherits(result, "plotly"))
})

