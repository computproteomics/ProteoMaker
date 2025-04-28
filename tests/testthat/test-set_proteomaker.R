library(testthat)
test_that("set_proteomaker returns correct default configuration", {
    config <- set_proteomaker()

    expect_equal(grepl("Proteome", config$fastaFilePath), TRUE)
    expect_equal(config$resultFilePath, "SimulatedDatasets")
    expect_equal(config$cores, 2)
    expect_equal(config$clusterType, "PSOCK")
    expect_equal(config$runStatTests, TRUE)
    expect_equal(config$calcAllBenchmarks, TRUE)
})

test_that("set_proteomaker allows custom configuration", {
    config <- set_proteomaker(fastaFilePath = "CustomProteomes", resultFilePath = "Results",
                           cores = 4, clusterType = "FORK", runStatTests = FALSE,
                           calcAllBenchmarks = FALSE)

    expect_equal(config$fastaFilePath, "CustomProteomes")
    expect_equal(config$resultFilePath, "Results")
    expect_equal(config$cores, 4)
    expect_equal(config$clusterType, "FORK")
    expect_equal(config$runStatTests, FALSE)
    expect_equal(config$calcAllBenchmarks, FALSE)
})

test_that("set_proteomaker creates the result directory", {
    temp_dir <- tempdir()
    config <- set_proteomaker(resultFilePath = file.path(temp_dir, "TestResults"))

    expect_true(dir.exists(config$resultFilePath))

    # Cleanup
    unlink(config$resultFilePath, recursive = TRUE)
})
