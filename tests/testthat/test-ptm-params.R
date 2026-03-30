library(testthat)

# Helper: build the pars list expected by modify_seq
make_pars <- function(lambda = 0) {
  pmod_res <- list(ph = c("S", "T", "Y"))
  ptms <- "ph"
  # pmod_res_distr is already multiplied by ptms_distr (as done inside modify())
  pmod_res_distr <- list(ph = c(S = 0.86, T = 0.13, Y = 0.01))
  ptms_distr <- list(ph = 1)
  list(pmod_res, ptms, pmod_res_distr, ptms_distr, lambda)
}

# Helper: minimal PTM parameters for run_sims
make_ptm_params <- function(params, lambda = 0) {
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

test_that("modify_seq with lambda = 0 yields approximately 1 modification per PTM type", {
  set.seed(42)
  # A sequence with many S, T, Y residues
  seq <- "MSTYSTSYTSSTTYYSSTMSTYSTSYTSSTTYYSS"
  pars <- make_pars(lambda = 0)

  results <- replicate(20, {
    r <- ProteoMaker:::modify_seq(seq, pars)
    length(r$Positions)
  })

  # lambda = 0 is clamped to 1e-4 internally; the truncated Poisson with such a
  # tiny lambda always returns 1, so each call should produce >= 1 modification.
  expect_true(all(results >= 1), info = paste("Modification counts:", paste(results, collapse = ", ")))
})

test_that("modify_seq with lambda = 0.0001 behaves the same as lambda = 0", {
  set.seed(42)
  seq <- "MSTYSTSYTSSTTYYSSTMSTYSTSYTSSTTYYSS"
  pars_0 <- make_pars(lambda = 0)
  pars_tiny <- make_pars(lambda = 0.0001)

  # Both should yield the same distribution of modification counts since
  # 0.0001 == 1e-4 and the internal clamp treats them identically.
  set.seed(7); counts_0 <- replicate(20, length(ProteoMaker:::modify_seq(seq, pars_0)$Positions))
  set.seed(7); counts_tiny <- replicate(20, length(ProteoMaker:::modify_seq(seq, pars_tiny)$Positions))

  expect_identical(counts_0, counts_tiny)
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

test_that("run_sims with PTMs and lambda=0 completes without error", {
  ll <- list.files(tempdir(), pattern = "output", full.names = TRUE)
  unlink(ll, recursive = TRUE)

  config <- test_proteomaker_config(resultFilePath = tempdir())
  params <- def_param()
  params$paramGroundTruth$NumReps <- 2
  params$paramGroundTruth$NumCond <- 2
  params <- make_ptm_params(params, lambda = 0)

  results <- run_sims(params, config)

  expect_true(is.list(results))
  expect_true(length(results) > 0)
})

test_that("flat PTMTypesDistr format list(ph=1) runs without error", {
  ll <- list.files(tempdir(), pattern = "output", full.names = TRUE)
  unlink(ll, recursive = TRUE)

  config <- test_proteomaker_config(resultFilePath = tempdir())
  params <- def_param()
  params$paramGroundTruth$NumReps <- 2
  params$paramGroundTruth$NumCond <- 2
  # Use the flat format documented in the YAML description
  params$paramGroundTruth$FracModProt <- 0.5
  params$paramGroundTruth$PropModPerProt <- 1
  params$paramGroundTruth$RemoveNonModFormFrac <- 0
  params$paramGroundTruth$PTMTypes <- list(mods = c("ph"))
  params$paramGroundTruth$PTMTypesMass <- list(ph = 79.966331)
  params$paramGroundTruth$PTMTypesDistr <- list(ph = 1)
  params$paramGroundTruth$ModifiableResidues <- list(mods = list(ph = c("S", "T", "Y")))
  params$paramGroundTruth$ModifiableResiduesDistr <- list(mods = list(ph = c(S = 0.86, T = 0.13, Y = 0.01)))
  params$paramGroundTruth$PTMMultipleLambda <- 0

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
  params <- make_ptm_params(params, lambda = 0)

  results <- run_sims(params, config)

  gt <- get_simulation(results[[1]]$Param, config, stage = "GroundTruth")
  gt_df <- gt$groundTruth

  # Every protein should appear exactly twice: 1 modified + 1 unmodified proteoform
  counts_per_protein <- table(gt_df$Accession)
  expect_true(all(counts_per_protein == 2),
    info = paste("Proteoform counts per protein:",
      paste(unique(as.integer(counts_per_protein)), collapse = ", ")))
})
