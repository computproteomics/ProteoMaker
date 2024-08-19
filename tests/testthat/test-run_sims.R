library(testthat)

test_that("run_sims completes without errors for minimal input", {
    params <- def_param()
    ll <- list.files(tempdir(), pattern="output", full.names = TRUE)
    unlink(ll, recursive = TRUE)
    config <- set_phosfake(resultFilePath = tempdir())
    params$paramGroundTruth$NumReps <- 2
    params$paramGroundTruth$NumCond <- 2
    params$paramDigest$DigestionEnzyme <- "Trypsin"

    results <- run_sims(params, config)
    
    expect_true(is.list(results))
    expect_true(length(results) > 0)
})

test_that("run_sims creates expected result files", {
    params <- def_param()
    ll <- list.files(tempdir(), pattern="output", full.names = TRUE)
    unlink(ll, recursive = TRUE)
    params$paramGroundTruth$NumReps <- 2
    params$paramGroundTruth$NumCond <- 2
    params$paramDigest$DigestionEnzyme <- "Trypsin"

    ll <- list.files(tempdir(), pattern="output", full.names = TRUE)
    unlink(ll, recursive = TRUE)
    
    config <- set_phosfake(resultFilePath = tempdir())
    
    results <- run_sims(params, config)
    
    expected_files <- c("outputGroundTruth_", "outputProteoformAb_", "outputDigest_", "outputMSRun_", "outputDataAnalysis_")
    files_created <- list.files(config$resultFilePath)
    
    for (file_prefix in expected_files) {
        expect_true(any(grepl(file_prefix, files_created)))
    }
})
