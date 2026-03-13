test_proteomaker_config <- function(...) {
  dots <- list(...)
  has_user_cores <- "cores" %in% names(dots)
  if (!has_user_cores) {
    dots$cores <- 2
  }
  if (is.null(dots$clusterType)) {
    dots$clusterType <- if (.Platform$OS.type == "unix") "FORK" else "PSOCK"
  }
  config <- do.call(set_proteomaker, dots)

  # In restricted R CMD check sandboxes, socket-based clusters may be blocked.
  # Keep parallel defaults for normal test runs, but force serial in check unless explicitly overridden.
  in_r_cmd_check <- nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))
  force_parallel <- identical(Sys.getenv("PROTEOMAKER_TEST_FORCE_PARALLEL"), "1")
  if (in_r_cmd_check && !force_parallel && !has_user_cores) {
    config$cores <- NULL
  }

  config
}
