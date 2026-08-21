# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(ProteoMaker)

cleanup_outputs <- function() {
  dirs <- c("Results", "SimulatedDatasets")
  unlink(c(dirs, file.path("..", dirs)), recursive = TRUE)
}
cleanup_outputs()
on.exit(cleanup_outputs(), add = TRUE)

test_check("ProteoMaker")
