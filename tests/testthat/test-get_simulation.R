library(testthat)

test_that("get_simulation returns NULL if the file does not exist", {
    config <- set_phosfake(resultFilePath = tempdir())
    param <- list(NumReps = 2, NumCond = 3)
    
    result <- get_simulation(param, config, stage = "MSRun")
    expect_null(result)
})

test_that("get_simulation returns the correct loaded data", {
    config <- set_phosfake(resultFilePath = tempdir())
    Param <- def_param()
    Param$paramGroundTruth$NumReps <- 5
    benchmarks <- run_sims(Param, config)

    result <- get_simulation(benchmarks[[1]]$Param, config, stage="MSRun")
    expect_true(is.list(result))
    expect_named(result, c("AfterMSRun", "Param"))
    
    result <- get_simulation(benchmarks[[1]]$Param, config)
    expect_true(is.list(result))
    expect_named(result, c("Benchmarks", "Param", "Stats", "StatsPep"))
})

test_that("get_simulation handles invalid stage name", {
    config <- set_phosfake(resultFilePath = tempdir())
    param <- list(NumReps = 2, NumCond = 3)
    
    expect_error(get_simulation(param, config, stage = "InvalidStage"), "Invalid stage name")
})
