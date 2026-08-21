## Regression test: setting an attribute on NULL when a protein group has fewer
## peptidoforms than MinUniquePep used to crash parLapply with
## "attempt to set an attribute on NULL".

make_peptable_below_threshold <- function() {
  # P1 has 3 peptidoforms (>= MinUniquePep = 3) -> summarized
  # P2 has 1 peptidoform  (<  MinUniquePep = 3) -> NULL from summarizeProtein
  data.frame(
    Accession = I(list(c("P1"), c("P1"), c("P1"), c("P2"))),
    PTMType   = I(replicate(4L, character(0), simplify = FALSE)),
    Sequence  = c("SEQ1", "SEQ2", "SEQ3", "SEQ4"),
    C_1_R_1   = c(1, 2, 3, 4),
    C_1_R_2   = c(2, 3, 4, 5),
    stringsAsFactors = FALSE
  )
}

test_that("proteinSummarisation does not error when a group is below MinUniquePep (sequential)", {
  peptable <- make_peptable_below_threshold()
  params <- list(
    ProtSummarization = "mean",
    MinUniquePep      = 3,
    IncludeModPep     = TRUE,
    SharedPep         = TRUE,
    QuantColnames     = c("C_1_R_1", "C_1_R_2"),
    Cores             = NULL
  )
  expect_no_error(res <- ProteoMaker:::proteinSummarisation(peptable, params))
  # Only P1 survives the threshold filter
  expect_equal(nrow(res), 1L)
  expect_equal(rownames(res), "P1")
})

test_that("proteinSummarisation does not error when a group is below MinUniquePep (parallel PSOCK)", {
  skip_on_cran()
  socket <- try(parallel::serverSocket(0L), silent = TRUE)
  skip_if(inherits(socket, "try-error"), "parallel sockets unavailable")
  close(socket)
  peptable <- make_peptable_below_threshold()
  params <- list(
    ProtSummarization = "mean",
    MinUniquePep      = 3,
    IncludeModPep     = TRUE,
    SharedPep         = TRUE,
    QuantColnames     = c("C_1_R_1", "C_1_R_2"),
    Cores             = 2L,
    ClusterType       = "PSOCK"
  )
  expect_no_error(res <- ProteoMaker:::proteinSummarisation(peptable, params))
  expect_equal(nrow(res), 1L)
  expect_equal(rownames(res), "P1")
})
