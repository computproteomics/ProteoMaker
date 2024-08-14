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