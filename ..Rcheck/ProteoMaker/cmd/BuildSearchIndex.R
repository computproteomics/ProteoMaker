#!/usr/bin/env Rscript

# Build a peptide search index that is identical to the one used in the
# ProteoMaker pipeline. Instead of reimplementing the logic we load the
# package, read the parameter YAML, and call the same helper that the
# digestion stage invokes (buildSearchIndexFromFasta).

suppressPackageStartupMessages({
  library(ProteoMaker)
})

args <- commandArgs(trailingOnly = TRUE)

yaml_path <- if (length(args) >= 1) args[[1]] else
  system.file("config", "parameters.yaml", package = "ProteoMaker")
out_rds   <- if (length(args) >= 2) args[[2]] else NA_character_
fasta_base <- if (length(args) >= 3) args[[3]] else
  system.file("Proteomes", package = "ProteoMaker")

if (!file.exists(yaml_path)) {
  stop("Parameter YAML not found: ", yaml_path)
}

# Use the package helper so that parameter handling stays in sync with the
# pipeline (defaults, choices, etc.).
param_categories <- ProteoMaker::def_param(yaml_path)

flatten_params <- function(cats) {
  parts <- c(
    cats$paramGroundTruth,
    cats$paramProteoformAb,
    cats$paramDigest,
    cats$paramMSRun,
    cats$paramDataAnalysis
  )
  # Remove NULL entries that arise when a category is missing
  parts[!vapply(parts, is.null, logical(1))]
}

idx_param <- flatten_params(param_categories)

if (!"PathToFasta" %in% names(idx_param) || is.na(idx_param$PathToFasta) ||
      !nzchar(idx_param$PathToFasta)) {
  stop("PathToFasta is not defined in the provided parameters.")
}

is_absolute_path <- function(path) {
  grepl("^(?:[A-Za-z]:|/|\\\\\\\\)", path)
}

fasta_path <- idx_param$PathToFasta
if (!is_absolute_path(fasta_path)) {
  fasta_path <- file.path(fasta_base, fasta_path)
}

if (!file.exists(fasta_path)) {
  stop("FASTA file not found: ", fasta_path)
}

idx_param$PathToFasta <- fasta_path

# Ensure cluster/cores fields exist (the helper expects them in some contexts)
if (!"Cores" %in% names(idx_param)) idx_param$Cores <- 1L
if (!"ClusterType" %in% names(idx_param)) idx_param$ClusterType <- "PSOCK"

cat("\n#PROTEOME INDEX - Using parameters from:", yaml_path, "\n")
cat(" + FASTA:", idx_param$PathToFasta, "\n")

search_index <- ProteoMaker:::buildSearchIndexFromFasta(idx_param)

if (!is.na(out_rds) && nzchar(out_rds)) {
  dir.create(dirname(out_rds), showWarnings = FALSE, recursive = TRUE)
  saveRDS(search_index, file = out_rds)
  cat("\nSaved index to:", out_rds, "\n")
}

cat("\n#PROTEOME INDEX - Done\n")
