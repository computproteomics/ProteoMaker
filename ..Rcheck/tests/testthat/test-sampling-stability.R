test_that("Donor sampling stable under extreme/absent log counts", {
  proteins <- data.frame(
    Accession = c("PX1"),
    Sequence = c("MKRPEPTIDERKPEPTK"),
    stringsAsFactors = FALSE
  )
  params <- list(
    Enzyme = "trypsin",
    PepMinLength = 6,
    PepMaxLength = 12,
    MaxNumMissedCleavages = 1,
    PTMTypes = list(c("ph", "ox")),
    ModifiableResidues = list(list(ph = c("S", "T", "Y"), ox = c("M")))
  )
  idx <- ProteoMaker:::buildSearchIndexFromSequences(proteins, params)
  e1 <- idx$proteins[[1]]
  L <- length(e1$win_start)
  skip_if(L < 2, "not enough windows in tiny index")

  # Extreme, equal log-counts -> uniform, should not overflow/underflow
  idx$proteins[[1]]$win_log_count <- rep(1000, L)
  donors <- ProteoMaker:::.pm_sample_uniform_peptidoforms(idx, params, N = min(10, L))
  expect_s3_class(donors, "data.frame")
  expect_equal(nrow(donors), min(10, L))

  # Missing log counts -> fallback to wc in log-space; still stable
  idx$proteins[[1]]$win_log_count <- NULL
  donors2 <- ProteoMaker:::.pm_sample_uniform_peptidoforms(idx, params, N = min(7, L))
  expect_s3_class(donors2, "data.frame")
  expect_equal(nrow(donors2), min(7, L))
})
