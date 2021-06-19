#' Composite likelihood for estimating LIR(Lagged Identification Rate)
#'
#' This CL function requires a sample format of pairwise lagged identification,
#' which is a list-like data containing lagged identification of of all pairs of
#' observation.
#'
#' This function is recommended when total observation time is very small and
#' pairwise lagged identification is easy to pre-calculate through
#' @seealso [move.pairwise()]. The notation of the estimation is \eqn{\hat{R_{tau}}}
#'
#' @param theta Parameter to calculate R_tau
#' @param model Model to calculate R_tau
#' @param n Number of individuals at each observation
#' @param m List of lagged identification for all observation pair
#' @param ... extra arguments to be passed to model to calculate R_tau
#'
#' @return Likelihood score(numeric). Inf if theta out of [0,1]
#' @export
#'
#' @examples
move.CL.pair <- function(theta, model, n, m, ...) {
  R_tau = model(theta, ...)
  if (min(R_tau * n)<=0) return(Inf);
  if (max(R_tau * n)>=1) return(Inf);
  object <- sum(m * log(R_tau*n) + (n-m) * log(1 - R_tau*n))
  return(-object)
}

move.MCLE <-
  function(theta,
           model,
           n,
           N,
           data,
           ...,
           optimizer = NULL,
           opt_arg = NULL) {

  }

#' LIR model A. The population is constant and no migrations occur.
#'
#' @param theta here theta is the estimated LIR(R_tau) directly
#'
#' @return theta
#' @export
#'
#' @examples
move.model.A <- function(theta) {
  return(theta)
}

#' LIR model B.
#'
#' @param theta (\eqn{\alpha, \beta})
#' @param tau List of time interval of lagged identification
#'
#' @return List of estimated R_tau
#' @export
#'
#' @examples
move.model.B <- function(theta, tau) {
  return(theta[1] * exp(-theta[2]*tau))
}

#' LIR model C.
#'
#' @param theta
#' @param tau
#'
#' @return
#' @export
#'
#' @examples
move.model.C <- function(theta, tau) {
  return(theta[1]*exp(-theta[2]*tau) + theta[3])
}
