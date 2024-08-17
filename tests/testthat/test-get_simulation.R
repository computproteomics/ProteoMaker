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
    Param$paramGroundTruth$NumReps <- c(3)
    benchmarks <- run_sims(Param, config)

    result <- get_simulation(benchmarks[[1]]$Param, config, stage="MSRun")
    expect_true(is.list(result))
    expect_named(result, c("AfterMSRun", "Param"))
    result <- get_simulation(benchmarks[[1]]$Param, config)
    expect_true(is.list(result))
    expect_named(result, c("Benchmarks", "Param", "Stats", "StatsPep"))
    result <- get_simulation(benchmarks[[1]]$Param, config, stage="GroundTruth")
    expect_true(is.list(result))
    expect_named(result, c("groundTruth", "Param"))
    result <- get_simulation(benchmarks[[1]]$Param, config, stage="Digest")
    expect_true(is.list(result))
    expect_named(result, c("BeforeMS", "Param"))
    result <- get_simulation(benchmarks[[1]]$Param, config, stage="ProteoformAb")
    expect_true(is.list(result))
    expect_named(result, c("Param", "proteoformAb"))
    # delete files in tempdir
    unlink(tempdir(), recursive = TRUE)
})

test_that("get_simulation returns the correct loaded data, no stats", {
    config <- set_phosfake(resultFilePath = tempdir(), runStatTests = FALSE)
    Param <- def_param()
    Param$paramGroundTruth$NumReps <- 5
    benchmarks <- run_sims(Param, config)
    
    result <- get_simulation(benchmarks[[1]]$Param, config, stage="MSRun")
    expect_true(is.list(result))
    expect_named(result, c("AfterMSRun", "Param"))
    result <- get_simulation(benchmarks[[1]]$Param, config)
    expect_true(is.null(result))
    result <- get_simulation(benchmarks[[1]]$Param, config, stage="GroundTruth")
    expect_true(is.list(result))
    expect_named(result, c("groundTruth", "Param"))
    result <- get_simulation(benchmarks[[1]]$Param, config, stage="Digest")
    expect_true(is.list(result))
    expect_named(result, c("BeforeMS", "Param"))
    result <- get_simulation(benchmarks[[1]]$Param, config, stage="ProteoformAb")
    expect_true(is.list(result))
    expect_named(result, c("Param", "proteoformAb"))
    # delete files in tempdir
    unlink(tempdir(), recursive = TRUE)
})

test_that("get_simulation returns the correct loaded data, no benchmarks", {
    config <- set_phosfake(resultFilePath = tempdir(), calcAllBenchmarks = FALSE)
    Param <- def_param()
    Param$paramGroundTruth$NumReps <- 5
    benchmarks <- run_sims(Param, config)
    
    result <- get_simulation(benchmarks[[1]]$Param, config, stage="MSRun")
    expect_true(is.list(result))
    expect_named(result, c("AfterMSRun", "Param"))
    result <- get_simulation(benchmarks[[1]]$Param, config)
    expect_true(is.list(result))
    expect_named(result, c("Benchmarks", "Param", "Stats", "StatsPep"))
    expect_true(is.null(result$Benchmarks))
    result <- get_simulation(benchmarks[[1]]$Param, config, stage="GroundTruth")
    expect_true(is.list(result))
    expect_named(result, c("groundTruth", "Param"))
    result <- get_simulation(benchmarks[[1]]$Param, config, stage="Digest")
    expect_true(is.list(result))
    expect_named(result, c("BeforeMS", "Param"))
    result <- get_simulation(benchmarks[[1]]$Param, config, stage="ProteoformAb")
    expect_true(is.list(result))
    expect_named(result, c("Param", "proteoformAb"))
    # delete files in tempdir
    unlink(tempdir(), recursive = TRUE)
})

test_that("get_simulation returns the correct loaded data, no stats, not benchmarks", {
    config <- set_phosfake(resultFilePath = tempdir(), calcAllBenchmarks = FALSE, runStatTests = FALSE)
    Param <- def_param()
    Param$paramGroundTruth$NumReps <- 5
    benchmarks <- run_sims(Param, config)
    
    result <- get_simulation(benchmarks[[1]]$Param, config, stage="MSRun")
    expect_true(is.list(result))
    expect_named(result, c("AfterMSRun", "Param"))
    result <- get_simulation(benchmarks[[1]]$Param, config)
    expect_true(is.null(result))
    result <- get_simulation(benchmarks[[1]]$Param, config, stage="GroundTruth")
    expect_true(is.list(result))
    expect_named(result, c("groundTruth", "Param"))
    result <- get_simulation(benchmarks[[1]]$Param, config, stage="Digest")
    expect_true(is.list(result))
    expect_named(result, c("BeforeMS", "Param"))
    result <- get_simulation(benchmarks[[1]]$Param, config, stage="ProteoformAb")
    expect_true(is.list(result))
    expect_named(result, c("Param", "proteoformAb"))
    # delete files in tempdir
    unlink(tempdir(), recursive = TRUE)
})




test_that("get_simulation handles invalid stage name", {
    config <- set_phosfake(resultFilePath = tempdir())
    param <- list(NumReps = 2, NumCond = 3)
    
    expect_error(get_simulation(param, config, stage = "InvalidStage"), "Invalid stage name")
})
