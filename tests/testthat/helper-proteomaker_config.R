test_proteomaker_config <- function(...) {
  cluster_type <- if (.Platform$OS.type == "unix") "FORK" else "PSOCK"
  set_proteomaker(..., clusterType = cluster_type)
}
