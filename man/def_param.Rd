% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/00_BatchRunFuncs.R
\name{def_param}
\alias{def_param}
\title{Generate default parameters or read from yaml file}
\usage{
def_param(yaml_file = NULL)
}
\arguments{
\item{yaml_file}{A character string specifying the path to the YAML file.
NULL means setting default values.}
}
\value{
A list containing categorized default parameters:
\describe{
  \item{paramGroundTruth}{Parameters related to ground truth data generation.}
  \item{paramProteoformAb}{Parameters related to proteoform abundance.}
  \item{paramDigest}{Parameters related to enzymatic digestion.}
  \item{paramMSRun}{Parameters related to the mass spectrometry run.}
  \item{paramDataAnalysis}{Parameters related to data analysis.}
}
}
\description{
This function reads a YAML file containing parameter settings and extracts
them into a structured list categorized by their types. The YAML file should
have a specific structure, with parameters organized under a top-level
`params` key, and each parameter containing metadata including its type
and default value.
}
\examples{
# Read YAML file from inst folder
yaml_path <- system.file("config", "params.yaml", package = "ProteoMaker")
params <- def_param(yaml_path)
}
