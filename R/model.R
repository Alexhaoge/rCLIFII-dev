#' @name LIR.model.A/B/C
#' @title LIR model A/B/C
#'
#' @description
#' A: \eqn{R_{\tau}=\alpha, \theta=\alpha=1/N}
#'
#' B: \eqn{R_{\tau}=\alpha e^{-\beta\tau}, \theta=c(\alpha, \beta)=c(1/N, \lambda)}
#'
#' C: \eqn{R_{\tau}=\gamma e^{-\beta\tau}+\alpha, \theta=c(\gamma, \beta, \alpha)=c(\frac{\lambda}{(\lambda+\mu)N}, \lambda+\mu, \frac{\mu}{(\lambda+\mu)N})}
#'
#' @param theta Parameter to estimate.
#' @param tau Lagged time \eqn{\tau}.
#'
#' @details
#' In model A, it is assumed that no migration occurs and population N remains constant.
#'
#' In model B, move-in rate equals to move-out rate so the population remains constant.
#' Note that migration in this model is permanent, so animals previously move out(or dead)
#' will not return.
#'
#' In model C, animals in the study area move out with probability of \eqn{\lambda} per
#' unit time and move in with probability of \eqn{\mu} per unit time. The population of the
#' whole area T is assumed to be constant. If \eqn{\lambda = \frac{\mu (Z-N)}{N}}, the
#' population within the study has an expectation of N, otherwise a warning will
#' be raised.
#'
#' @return \eqn{R_{\tau}}
#' @rdname model
#' @export
#'
#' @examples
LIR.model.A <- function(theta, tau, ...) {
  return(theta)
}

#' @name LIR.grad.A/B/C
#' @title Gradient of model A/B/C
#'
#' @param theta Parameter to estimate.
#' @param tau Lagged time \eqn{\tau}.
#'
#' @return vector
#' @export
#' @rdname grad
#'
#' @examples
LIR.grad.A <- function(theta, tau, ...) {
  return(c(1))
}

#' @name LIR.hessian.A/B/C
#' @title Hessian matrix of model A/B/C
#'
#' @param theta Parameter to estimate.
#' @param tau Lagged time \eqn{\tau}.
#'
#' @return A symmetric matrix.
#' @rdname hessian
#' @export
#'
#' @examples
LIR.hessian.A <- function(theta, tau, ...) {
  return(matrix(0, nrow = 1, ncol = 1))
}

#' @rdname model
#' @export
LIR.model.B <- function(theta, tau, ...) {
  return(theta[1] * exp(-theta[2]*tau))
}

#' @rdname grad
#' @export
LIR.grad.B <- function(theta, tau, ...) {
  return(c(1, -theta[1] * tau) * exp(-theta[2]*tau))
}

#' @rdname hessian
#' @export
LIR.hessian.B <- function(theta, tau, ...) {
  return(matrix(
    c(0, -tau, -tau, theta[1] * tau^2) * exp(-theta[2] * tau),
    nrow = 2,
    ncol = 2
  ))
}

#' @rdname model
#' @export
LIR.model.C <- function(theta, tau, ...) {
  return(theta[1]*exp(-theta[2]*tau) + theta[3])
}

#' @rdname grad
#' @export
LIR.grad.C <- function(theta, tau, ...) {
  return(c(exp(-theta[2]*tau), -theta[1] * tau * exp(-theta[2]*tau), 1))
}

#' @rdname hessian
LIR.hessian.C <- function(theta, tau, ...) {
  return(matrix(
    c(0, -tau, 0, -tau, theta[1] * tau^2, 0, 0, 0, 0) * exp(-theta[2] * tau),
    nrow = 3,
    ncol = 3
  ))
}
