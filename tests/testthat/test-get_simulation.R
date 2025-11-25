library(testthat)

test_that("get_simulation returns NULL if the file does not exist", {
  config <- set_proteomaker(resultFilePath = tempdir())
  param <- list(NumReps = 2, NumCond = 3)

  result <- get_simulation(param, config, stage = "MSRun")
  expect_null(result)
})

test_that("get_simulation returns the correct loaded data", {
  ll <- list.files(tempdir(), pattern = "output", full.names = TRUE)
  unlink(ll, recursive = TRUE)
  config <- set_proteomaker(resultFilePath = tempdir())
  Param <- def_param()
  Param$paramGroundTruth$NumReps <- c(3)
  benchmarks <- run_sims(Param, config)

  result <- get_simulation(benchmarks[[1]]$Param, config, stage = "MSRun")
  expect_true(is.list(result))
  expect_named(result, c("AfterMSRun", "Param"), ignore.order = TRUE)
  result <- get_simulation(benchmarks[[1]]$Param, config)
  expect_true(is.list(result))
  expect_named(result, c("Benchmarks", "Param", "Stats", "StatsPep"), ignore.order = TRUE)
  result <- get_simulation(benchmarks[[1]]$Param, config, stage = "GroundTruth")
  expect_true(is.list(result))
  expect_named(result, c("groundTruth", "Param"), ignore.order = TRUE)
  result <- get_simulation(benchmarks[[1]]$Param, config, stage = "Digest")
  expect_true(is.list(result))
  expect_named(result, c("BeforeMS", "Param"), ignore.order = TRUE)
  result <- get_simulation(benchmarks[[1]]$Param, config, stage = "ProteoformAb")
  expect_true(is.list(result))
  expect_named(result, c("Param", "proteoformAb"), ignore.order = TRUE)
})

test_that("get_simulation returns the correct loaded data, no stats", {
  ll <- list.files(tempdir(), pattern = "output", full.names = TRUE)
  unlink(ll, recursive = TRUE)

  config <- set_proteomaker(resultFilePath = tempdir(), runStatTests = FALSE)
  Param <- def_param()
  Param$paramGroundTruth$NumReps <- 5
  benchmarks <- run_sims(Param, config)

  result <- get_simulation(benchmarks[[1]]$Param, config, stage = "MSRun")
  expect_true(is.list(result))
  expect_named(result, c("AfterMSRun", "Param"), ignore.order = TRUE)
  result <- get_simulation(benchmarks[[1]]$Param, config)
  expect_true(is.null(result))
  result <- get_simulation(benchmarks[[1]]$Param, config, stage = "GroundTruth")
  expect_true(is.list(result))
  expect_named(result, c("groundTruth", "Param"), ignore.order = TRUE)
  result <- get_simulation(benchmarks[[1]]$Param, config, stage = "Digest")
  expect_true(is.list(result))
  expect_named(result, c("BeforeMS", "Param"), ignore.order = TRUE)
  result <- get_simulation(benchmarks[[1]]$Param, config, stage = "ProteoformAb")
  expect_true(is.list(result))
  expect_named(result, c("Param", "proteoformAb"), ignore.order = TRUE)
})

test_that("get_simulation returns the correct loaded data, no benchmarks", {
  ll <- list.files(tempdir(), pattern = "output", full.names = TRUE)
  unlink(ll, recursive = TRUE)

  config <- set_proteomaker(resultFilePath = tempdir(), calcAllBenchmarks = FALSE)
  Param <- def_param()
  Param$paramGroundTruth$NumReps <- 5
  benchmarks <- run_sims(Param, config)

  result <- get_simulation(benchmarks[[1]]$Param, config, stage = "MSRun")
  expect_true(is.list(result))
  expect_named(result, c("AfterMSRun", "Param"), ignore.order = TRUE)
  result <- get_simulation(benchmarks[[1]]$Param, config)
  expect_true(is.list(result))
  expect_named(result, c("Benchmarks", "Param", "Stats", "StatsPep"), ignore.order = TRUE)
  expect_true(is.null(result$Benchmarks))
  result <- get_simulation(benchmarks[[1]]$Param, config, stage = "GroundTruth")
  expect_true(is.list(result))
  expect_named(result, c("groundTruth", "Param"), ignore.order = TRUE)
  result <- get_simulation(benchmarks[[1]]$Param, config, stage = "Digest")
  expect_true(is.list(result))
  expect_named(result, c("BeforeMS", "Param"), ignore.order = TRUE)
  result <- get_simulation(benchmarks[[1]]$Param, config, stage = "ProteoformAb")
  expect_true(is.list(result))
  expect_named(result, c("Param", "proteoformAb"), ignore.order = TRUE)
})

test_that("get_simulation returns the correct loaded data, no stats, not benchmarks", {
  ll <- list.files(tempdir(), pattern = "output", full.names = TRUE)
  unlink(ll, recursive = TRUE)

  config <- set_proteomaker(resultFilePath = tempdir(), calcAllBenchmarks = FALSE, runStatTests = FALSE)
  Param <- def_param()
  Param$paramGroundTruth$NumReps <- 5
  benchmarks <- run_sims(Param, config)

  result <- get_simulation(benchmarks[[1]]$Param, config, stage = "MSRun")
  expect_true(is.list(result))
  expect_named(result, c("AfterMSRun", "Param"), ignore.order = TRUE)
  result <- get_simulation(benchmarks[[1]]$Param, config)
  expect_true(is.null(result))
  result <- get_simulation(benchmarks[[1]]$Param, config, stage = "GroundTruth")
  expect_true(is.list(result))
  expect_named(result, c("groundTruth", "Param"), ignore.order = TRUE)
  result <- get_simulation(benchmarks[[1]]$Param, config, stage = "Digest")
  expect_true(is.list(result))
  expect_named(result, c("BeforeMS", "Param"), ignore.order = TRUE)
  result <- get_simulation(benchmarks[[1]]$Param, config, stage = "ProteoformAb")
  expect_true(is.list(result))
  expect_named(result, c("Param", "proteoformAb"), ignore.order = TRUE)
})


test_that("get_simulation handles invalid stage name", {
  config <- set_proteomaker(resultFilePath = tempdir())
  param <- list(NumReps = 2, NumCond = 3)

  expect_error(get_simulation(param, config, stage = "InvalidStage"), "Invalid stage name")
})
