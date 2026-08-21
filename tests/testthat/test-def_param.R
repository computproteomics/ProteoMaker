library(testthat)

test_that("def_param loads parameters from YAML file", {
  # Use a temporary YAML file for testing
  temp_yaml <- tempfile(fileext = ".yaml")
  writeLines("
params:
  test_param:
    type: paramGroundTruth
    default: 42
  another_param:
    type: paramMSRun
    default: 3.14
", temp_yaml)

  params <- def_param(temp_yaml)

  expect_equal(params$paramGroundTruth$test_param, 42)
  expect_equal(params$paramMSRun$another_param, 3.14)

  # Cleanup
  unlink(temp_yaml)
})

test_that("def_param handles NA values correctly", {
  temp_yaml <- tempfile(fileext = ".yaml")
  writeLines("
params:
  test_param:
    type: paramGroundTruth
    default: 'NA'
  another_param:
    type: paramMSRun
    default: 3.14
", temp_yaml)

  params <- def_param(temp_yaml)

  expect_true(is.na(params$paramGroundTruth$test_param))
  expect_equal(params$paramMSRun$another_param, 3.14)

  # Cleanup
  unlink(temp_yaml)
})

test_that("def_param uses default YAML file if none provided", {
  # This test assumes you have a default YAML file in your package's inst/config directory.
  params <- def_param()

  # Check that params are non-null
  expect_true(length(params) > 0)
})

test_that("def_param reads configured PTM maps", {
  params <- def_param(system.file("config", "parameters_human_ph_ox_baseline.yaml",
                                  package = "ProteoMaker"))
  ptm <- params$paramGroundTruth

  expect_equal(ptm$PTMTypes$mods, c("ph", "ox"))
  expect_equal(names(ptm$PTMTypesDistr[[1]]), c("ph", "ox"))
  expect_equal(ptm$PTMTypesDistr[[1]]$ph, 0.5)
  expect_equal(ptm$PTMTypesMass[[1]]$ox, 15.994915)
  expect_equal(ptm$ModifiableResidues[[1]]$ph, c("S", "T", "Y"))
  expect_equal(ptm$ModifiableResiduesDistr[[1]]$ox, 1)
})
