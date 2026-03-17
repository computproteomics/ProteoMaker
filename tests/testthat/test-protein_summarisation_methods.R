make_peptable <- function(n_per_group = 3L) {
  stopifnot(n_per_group >= 1L)
  accessions <- c(rep("P1", n_per_group), rep("P2", n_per_group))
  seqs <- paste0("SEQ", seq_len(2L * n_per_group))
  data.frame(
    Accession = I(lapply(accessions, function(x) c(x))),
    PTMType = I(replicate(2L * n_per_group, character(0), simplify = FALSE)),
    Sequence = seqs,
    C_1_R_1 = seq_len(2L * n_per_group),
    C_1_R_2 = seq_len(2L * n_per_group) + 1,
    stringsAsFactors = FALSE
  )
}

make_peptable_missing <- function() {
  data.frame(
    Accession = I(list(c("P1"), c("P1"), c("P1"), c("P2"), c("P2"), c("P2"))),
    PTMType = I(replicate(6L, character(0), simplify = FALSE)),
    Sequence = paste0("MSEQ", seq_len(6L)),
    C_1_R_1 = c(1, NA, 3, 4, 5, NA),
    C_1_R_2 = c(2, 2, NA, NA, 6, 7),
    stringsAsFactors = FALSE
  )
}

test_that("proteinSummarisation supports common methods", {
  methods <- c("sum.top3", "median", "mean", "sum", "medpolish")
  peptable <- make_peptable(n_per_group = 3L)
  base_params <- list(
    MinUniquePep = 1,
    IncludeModPep = TRUE,
    SharedPep = TRUE,
    QuantColnames = c("C_1_R_1", "C_1_R_2"),
    Cores = NULL
  )

  for (m in methods) {
    params <- base_params
    params$ProtSummarization <- m
    res <- ProteoMaker:::proteinSummarisation(peptable, params)
    expect_equal(nrow(res), 2)
    expect_true(all(params$QuantColnames %in% colnames(res)))
    expect_true(all(!is.na(as.matrix(res[, params$QuantColnames]))))
  }
})

test_that("proteinSummarisation methods handle missing values", {
  methods <- c("sum.top3", "median", "mean", "sum", "medpolish")
  peptable <- make_peptable_missing()
  base_params <- list(
    MinUniquePep = 1,
    IncludeModPep = TRUE,
    SharedPep = TRUE,
    QuantColnames = c("C_1_R_1", "C_1_R_2"),
    Cores = NULL
  )

  for (m in methods) {
    params <- base_params
    params$ProtSummarization <- m
    res <- ProteoMaker:::proteinSummarisation(peptable, params)
    expect_gte(nrow(res), 1)
    expect_true(all(params$QuantColnames %in% colnames(res)))
    # At least one quantified value must remain after per-group NA filtering.
    expect_true(all(rowSums(is.na(res[, params$QuantColnames, drop = FALSE])) < length(params$QuantColnames)))
  }
})

test_that("proteinSummarisation medpolish works with >3 peptides per group", {
  peptable <- make_peptable(n_per_group = 4L)
  params <- list(
    ProtSummarization = "medpolish",
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
