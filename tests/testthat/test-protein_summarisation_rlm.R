test_that("proteinSummarisation rlm returns all samples including baseline", {
  peptable <- data.frame(
    Accession = I(list(c("P1"), c("P1"), c("P2"), c("P2"))),
    PTMType = I(list(character(0), character(0), character(0), character(0))),
    Sequence = c("AAA", "BBB", "CCC", "DDD"),
    C_1_R_1 = c(1, 2, 3, 4),
    C_1_R_2 = c(2, 3, 4, 5),
    stringsAsFactors = FALSE
  )

  params <- list(
    ProtSummarization = "rlm",
    MinUniquePep = 1,
    IncludeModPep = TRUE,
    SharedPep = TRUE,
    QuantColnames = c("C_1_R_1", "C_1_R_2"),
    Cores = NULL
  )

  res <- ProteoMaker:::proteinSummarisation(peptable, params)

  expect_equal(nrow(res), 2)
  expect_true(all(params$QuantColnames %in% colnames(res)))
  expect_true(all(!is.na(as.matrix(res[, params$QuantColnames]))))
})
