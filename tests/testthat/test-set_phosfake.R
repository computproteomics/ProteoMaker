library(testthat)
test_that("set_phosfake returns correct default configuration", {
    config <- set_phosfake()
    
    expect_equal(grepl("Proteome", config$fastaFilePath), TRUE)
    expect_equal(config$resultFilePath, "SimulatedDatasets")
    expect_equal(config$cores, 2)
    expect_equal(config$clusterType, "FORK")
    expect_equal(config$runStatTests, TRUE)
    expect_equal(config$calcAllBenchmarks, TRUE)
})

test_that("set_phosfake allows custom configuration", {
    config <- set_phosfake(fastaFilePath = "CustomProteomes", resultFilePath = "Results", 
                           cores = 4, clusterType = "PSOCK", runStatTests = FALSE, 
                           calcAllBenchmarks = FALSE)
    
    expect_equal(config$fastaFilePath, "CustomProteomes")
    expect_equal(config$resultFilePath, "Results")
    expect_equal(config$cores, 4)
    expect_equal(config$clusterType, "PSOCK")
    expect_equal(config$runStatTests, FALSE)
    expect_equal(config$calcAllBenchmarks, FALSE)
})

test_that("set_phosfake creates the result directory", {
    temp_dir <- tempdir()
    config <- set_phosfake(resultFilePath = file.path(temp_dir, "TestResults"))
    
    expect_true(dir.exists(config$resultFilePath))
    
    # Cleanup
    unlink(config$resultFilePath, recursive = TRUE)
})
