% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_Digestion.R
\name{addDetectability}
\alias{addDetectability}
\title{Calculate the detectability of a peptide sequence}
\usage{
addDetectability(peptides, parameters)
}
\arguments{
\item{peptides}{A character vector containing the peptide sequences.}

\item{parameters}{A list of parameters including the number of cores for parallel computing.}
}
\value{
A numeric vector representing the detectability of the peptides.
}
\description{
This function calculates the detectability of a peptide sequence based on the amino acid
using the PeptideRanger package. It returns a numeric vector representing the detectability
of the peptides. The prediction of the detectability is based on the amino acid composition
and does not take into account post-translational modifications.
}
\keyword{internal}
