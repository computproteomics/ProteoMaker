% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_Digestion.R
\name{filterDigestedProt}
\alias{filterDigestedProt}
\title{Create enriched and non-enriched fractions of proteolytic peptides}
\usage{
filterDigestedProt(DigestedProt, parameters)
}
\arguments{
\item{DigestedProt}{A data frame containing the digested proteolytic peptides and associated data.}

\item{parameters}{A list of parameters that includes EnrichmentLoss, ModificationLoss, EnrichmentEfficiency,
EnrichmentNoise, and QuantColnames.}
}
\value{
A list with two elements:
\describe{
  \item{NonEnriched}{A data frame containing the non-enriched peptide fraction, which includes both modified and non-modified peptides.}
  \item{Enriched}{A data frame containing the enriched peptide fraction, where modified peptides have been enriched based on the EnrichmentEfficiency, and noise has been added to simulate the enrichment process. If no modified peptides are present, this will be \code{NULL}.}
}
}
\description{
This function separates modified peptides into enriched and non-enriched fractions,
adjusts the peptide abundances based on enrichment efficiency, and introduces noise
due to the enrichment process. It returns a list containing the enriched and
non-enriched peptide sets.
}
