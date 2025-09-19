test_that("Index build and donor sampling are lightweight and consistent", {
  # Tiny protein set with simple tryptic sites
  proteins <- data.frame(
    Accession = c("A1", "A2"),
    Sequence = c("ACDEFGHIKLMKQR", "MPEPTIDERKPEPTK"),
    stringsAsFactors = FALSE
  )
  params <- list(
    Enzyme = "trypsin",
    PepMinLength = 7,
    PepMaxLength = 12,
    MaxNumMissedCleavages = 1,
    PTMTypes = list(c("ph", "ox")),
    ModifiableResidues = list(list(ph = c("S","T","Y"), ox = c("M")))
  )
  idx <- ProteoMaker:::buildSearchIndexFromSequences(proteins, params)
  expect_type(idx, "list")
  expect_true(length(idx$proteins) >= 1)
  # Each indexed protein should carry window vectors and log counts
  e1 <- idx$proteins[[1]]
  expect_true(length(e1$win_start) == length(e1$win_stop))
  expect_true(length(e1$win_start) == length(e1$win_log_count))
  # Sample a handful of donors and verify shape
  donors <- ProteoMaker:::.pm_sample_uniform_peptidoforms(idx, params, N = 8)
  expect_s3_class(donors, "data.frame")
  expect_equal(nrow(donors), 8)
  expect_true(all(c("Peptide","Start","Stop","MC") %in% names(donors)))
})

