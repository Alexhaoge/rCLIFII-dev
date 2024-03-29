% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{LIR.bootstrap}
\alias{LIR.bootstrap}
\title{Bootstrap on individuals}
\usage{
LIR.bootstrap(data, seed = NULL)
}
\arguments{
\item{data}{Observation matrix. Each row represent an individual
and each column represent an observation}

\item{seed}{random seed}
}
\value{
An new observation matrix. Same dims as the input data/
}
\description{
Perform a single bootstrap on the given observation matrix.
SRSWR will be performed on individuals (row of the matrix)
}
