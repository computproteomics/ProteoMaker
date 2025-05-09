% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/00_BatchRunFuncs.R
\name{visualize_one_sim}
\alias{visualize_one_sim}
\title{Visualize benchmarks for a specific simulation}
\usage{
visualize_one_sim(BenchMatrix, current_row = 1)
}
\arguments{
\item{BenchMatrix}{A data frame containing the benchmark results for multiple
simulations, typically generated by \code{matrix_benchmarks}.}

\item{current_row}{An integer or a character string specifying the row number
or the simulation identifier for which the benchmarks are to be visualized.
Default is 1.}
}
\value{
A plot object showing the visualized benchmarks and parameters for
the specified simulation.
}
\description{
This function creates visualizations of the benchmark results for a specific
simulation. It uses various plotting methods to display the normalized values
of benchmarks and parameters, aiding in the analysis of simulation outcomes.
}
\examples{
benchmarks <- matrix_benchmarks(run_sims(def_param(), set_proteomaker()), set_proteomaker())
visualize_one_sim(benchmarks, current_row = 1)
}
