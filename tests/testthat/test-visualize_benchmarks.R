library(testthat)

test_that("visualize_benchmarks returns a plotly object", {
    Param <- def_param()
    ll <- list.files(tempdir(), pattern="output", full.names = TRUE)
    unlink(ll, recursive = TRUE)
    proteomaker_config <- set_proteomaker(resultFilePath = tempdir())
    benchmarks <- run_sims(Param, proteomaker_config)
    benchmatrix <- matrix_benchmarks(benchmarks, proteomaker_config)
    result <- visualize_benchmarks(benchmatrix)
})

