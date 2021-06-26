#' LIR model A. The population is constant and no migrations occur.
#'
#' @param theta here theta is the estimated LIR(R_tau) directly
#' @param tau List of time interval of lagged identification
#'
#' @return rep(theta, length(tau))
#' @export
#'
#' @examples
LIR.model.A <- function(theta, tau) {
  return(theta)
}

LIR.grad.A <- function(theta, tau) {
  return(c(1))
}

LIR.hessian.A <- function(theta, tau) {
  return(matrix(0, nrow = 1, ncol = 1))
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
LIR.model.B <- function(theta, tau) {
  return(theta[1] * exp(-theta[2]*tau))
}

LIR.grad.B <- function(theta, tau) {
  return(c(1, -theta[1] * tau) * exp(-theta[2]*tau))
}

LIR.hessian.B <- function(theta, tau) {
  return(matrix(
    c(0, -tau, -tau, theta[1] * tau^2) * exp(-theta[2] * tau),
    nrow = 2,
    ncol = 2
  ))
}

#' LIR model C.
#'
#' @param theta (\eqn{\gamma, \beta, \alpha})
#' @param tau List of time interval of lagged identification
#'
#' @return
#' @export
#'
#' @examples
LIR.model.C <- function(theta, tau) {
  return(theta[1]*exp(-theta[2]*tau) + theta[3])
}

LIR.grad.C <- function(theta, tau) {
  return(c(exp(-theta[2]*tau), -theta[1] * tau * exp(-theta[2]*tau), 1))
}

LIR.hessian.C <- function(theta, tau) {
  return(matrix(
    c(0, -tau, 0, -tau, theta[1] * tau^2, 0, 0, 0, 0) * exp(-theta[2] * tau),
    nrow = 3,
    ncol = 3
  ))
}
