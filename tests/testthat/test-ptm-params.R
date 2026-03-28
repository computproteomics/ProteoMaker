library(testthat)

# Helper: build the pars list expected by modify_seq
make_pars <- function(lambda = 0.001) {
  pmod_res <- list(ph = c("S", "T", "Y"))
  ptms <- "ph"
  # pmod_res_distr is already multiplied by ptms_distr (as done inside modify())
  pmod_res_distr <- list(ph = c(S = 0.86, T = 0.13, Y = 0.01))
  ptms_distr <- list(ph = 1)
  list(pmod_res, ptms, pmod_res_distr, ptms_distr, lambda)
}

# Helper: minimal PTM parameters for run_sims
make_ptm_params <- function(params, lambda = 0.001) {
  params$paramGroundTruth$FracModProt <- 1
  params$paramGroundTruth$PropModPerProt <- 1
  params$paramGroundTruth$RemoveNonModFormFrac <- 0
  params$paramGroundTruth$PTMTypes <- list(mods = c("ph"))
  params$paramGroundTruth$PTMTypesMass <- list(ph = 79.966331)
  params$paramGroundTruth$PTMTypesDistr <- list(mods = c(ph = 1))
  params$paramGroundTruth$ModifiableResidues <- list(mods = list(ph = c("S", "T", "Y")))
  params$paramGroundTruth$ModifiableResiduesDistr <- list(mods = list(ph = c(S = 0.86, T = 0.13, Y = 0.01)))
  params$paramGroundTruth$PTMMultipleLambda <- lambda
  params
}

test_that("modify_seq with lambda = 0 runs without error for multi-residue sequence", {
  set.seed(42)
  seq <- "MSTYSTSYTSSTTYYSSTMSTYSTSYTSSTTYYSS"
  pars <- make_pars(lambda = 0)

  # Should not throw an error even when lambda = 0 and pos_len > 1,
  # and must produce no modifications (multi-site sampling is skipped).
  r <- ProteoMaker:::modify_seq(seq, pars)
  expect_equal(length(r$Positions), 0L)
})

test_that("modify_seq with small positive lambda produces at least one modification", {
  set.seed(42)
  # A sequence with many S, T, Y residues
  seq <- "MSTYSTSYTSSTTYYSSTMSTYSTSYTSSTTYYSS"
  pars <- make_pars(lambda = 0.001)

  results <- replicate(20, {
    r <- ProteoMaker:::modify_seq(seq, pars)
    length(r$Positions)
  })

  # With small positive lambda and many modifiable residues, modifications should occur
  expect_true(any(results >= 1), info = paste("Modification counts:", paste(results, collapse = ", ")))
})

test_that("modify_seq with lambda > 0 can produce multiple modifications", {
  set.seed(7)
  seq <- paste(rep("STYSTYSTYSTYSTYSTYSTYSTYS", 5), collapse = "")
  pars <- make_pars(lambda = 2)

  results <- replicate(20, {
    r <- ProteoMaker:::modify_seq(seq, pars)
    length(r$Positions)
  })

  # With a large lambda and many sites, we should occasionally see > 1 modification
  expect_true(any(results > 1), info = "Expected at least one run with > 1 modification")
})

test_that("run_sims with PTMs completes without error", {
  ll <- list.files(tempdir(), pattern = "output", full.names = TRUE)
  unlink(ll, recursive = TRUE)

  config <- test_proteomaker_config(resultFilePath = tempdir())
  params <- def_param()
  params$paramGroundTruth$NumReps <- 2
  params$paramGroundTruth$NumCond <- 2
  params <- make_ptm_params(params, lambda = 0.001)

  results <- run_sims(params, config)

  expect_true(is.list(results))
  expect_true(length(results) > 0)
})

test_that("PropModPerProt=1 and RemoveNonModFormFrac=0 yields 2 proteoforms per modified protein", {
  ll <- list.files(tempdir(), pattern = "output", full.names = TRUE)
  unlink(ll, recursive = TRUE)

  config <- test_proteomaker_config(resultFilePath = tempdir())
  params <- def_param()
  params$paramGroundTruth$NumReps <- 2
  params$paramGroundTruth$NumCond <- 2
  params <- make_ptm_params(params, lambda = 0.001)

  results <- run_sims(params, config)

  gt <- get_simulation(results[[1]]$Param, config, stage = "GroundTruth")
  gt_df <- gt$groundTruth

  # Every protein should appear exactly twice: 1 modified + 1 unmodified proteoform
  counts_per_protein <- table(gt_df$Accession)
  expect_true(all(counts_per_protein == 2),
    info = paste("Proteoform counts per protein:",
      paste(unique(as.integer(counts_per_protein)), collapse = ", ")))
})
