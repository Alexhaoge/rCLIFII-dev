% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{LIR.simulate.A/B/C}
\alias{LIR.simulate.A/B/C}
\alias{LIR.simulate.A}
\alias{LIR.simulate.B}
\alias{LIR.simulate.C}
\title{Simulate animal movement identification data with model A/B/C}
\usage{
LIR.simulate.A(N, n, tp, seed = NULL)

LIR.simulate.B(N, n, tp, lambda, seed = NULL)

LIR.simulate.C(Z, N, n, tp, lambda, mu, seed = NULL)
}
\arguments{
\item{N}{Population within the study area.(For model C, the initial population)}

\item{n}{Number of identification in each observation}

\item{tp}{Time of each observation if tp is a list or array, otherwise the number
of observation if t is an integer. Model B and C require tp to be a list.}

\item{seed}{random seed}

\item{lambda}{move-out rate}

\item{Z}{Total population of the whole area.(Only for model C)}

\item{mu}{move-in rate}
}
\value{
matrix of length(t) columns. The number of rows is N for model A/B,
Z for model C.
}
\description{
Simulate animal movement identification data with model A/B/C
}
\details{
In model A, it is assumed that no migration occurs and population N remains constant.

In model B, move-in rate equals to move-out rate so the population remains constant.
Note that migration in this model is permanent, so animals previously move out(or dead)
will not return.

In model C, animals in the study area move out with probability of \eqn{\lambda} per
unit time and move in with probability of \eqn{\mu} per unit time. The population of the
whole area T is assumed to be constant. If \eqn{\lambda = \frac{\mu (Z-N)}{N}}, the
population within the study has an expectation of N, otherwise a warning will
be raised.
}
\examples{
# Example of data simulation
# Set observation time
tp <- c(1:5, 51:55, 101:105, 501:505, 601:605)
# Generate observation matrix with model C
data <- LIR.simulate.C(Z=300, N=100, n=40, tp=tp, lambda=0.08, mu=0.04)
dim(data)
# [1] 300  25

}
