library(testthat)

test_that("matrix_benchmarks returns a correctly formatted data frame", {
    Param <- def_param()
    Param$paramGroundTruth$NumReps <- 5
    phosfake_config <- set_phosfake(resultFilePath = tempdir())
    benchmarks <- run_sims(Param, phosfake_config)
    results <- matrix_benchmarks(benchmarks, phosfake_config)

    expect_true(is.data.frame(results))
    expect_true(results$NumReps == 5)
    
    })
