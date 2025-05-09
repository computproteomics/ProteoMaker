% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/00_BatchRunFuncs.R
\name{plot_params}
\alias{plot_params}
\title{Plot parameters for a specific simulation}
\usage{
plot_params(BenchMatrix, current_row = 1)
}
\arguments{
\item{BenchMatrix}{A data frame containing the parameter values for multiple
simulations.}

\item{current_row}{An integer or a character string specifying the row number
or the simulation identifier for which the parameters are to be visualized.}
}
\value{
A plot object visualizing the parameter values for the specified
simulation.
}
\description{
This function generates visualizations of the parameter values for a specific
simulation. It uses a gauge-like representation to show the range and actual
values of parameters, providing an intuitive understanding of the simulation
setup.
}
\examples{
\dontrun{
benchmarks <- matrix_benchmarks(run_sims(def_param(), set_proteomaker()), set_proteomaker())
plot_params(benchmarks, current_row = 1)
}
}
