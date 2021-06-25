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
  R_tau <- model(theta, ...)
  if (min(R_tau * n)<=0) return(Inf);
  if (max(R_tau * n)>=1) return(Inf);
  object <- sum(m * log(R_tau*n) + (n-m) * log(1 - R_tau*n))
  return(-object)
}

#' TODO: undecided
move.CL <- function() {

}

#' Maximal Composite Likelihood Estimate for pairwise data
#'
#' @param theta Initial value for parameters to estimate
#' @param model Model to calculate R_tau
#' @param n Number of individuals at each observation
#' @param m List of lagged identification for all observation pair
#' @param tau List of time interval of lagged identification
#' @param ... Additional parameters passed to model
#' @param lower Lower bound for estimation
#' @param upper Upper bound for estimation
#' @param optimizer Optimizer to use. Default NULL(which means @seealso [optim()]
#'   will be used). If a user-defined optimizer is applied, it is highly recommended
#'   to wrap this optimizer so that it has the same parameter format as optim().
#' @param opt_arg Control list for optimizer. Default empty list
#' @param verbose Return detailed output. Default FALSE.
#'
#' @return If verbose is TRUE return the full result from optimizer otherwise return MCLE only.
#' @export
#'
#' @examples
move.MCLE.pair <-
  function(theta,
           model,
           n,
           m,
           tau,
           ...,
           lower = 0.0,
           upper = 1.0,
           optimizer = NULL,
           opt_arg = list(),
           verbose = FALSE) {
    clf <- function(theta_) {move.CL.pair(theta_, model = model, n = n, m = m, tau = tau, ...)}
    if (is.null(optimizer)) {
      opt_res <-
        optim(
          par = theta,
          fn = clf,
          method = ifelse(length(theta) == 1, 'Brent', 'SANN'),
          lower = lower,
          upper = upper,
          control = opt_arg
        )
    } else {
      opt_res <- optimizer(theta, clf, lower = lower, upper = upper, control = opt_arg)
    }
    return(ifelse(verbose, opt_res, opt_res$par))
  }

#' Maximal Composite Likelihood Estimate
#' If number of observation(length(t)) is less than 20000, pairwise MCLE
#' will be used. Otherwise pairwise lagged identification will be
#' calculated seperately at every iteration.
#'
#' @param theta Initial value for parameters to estimate
#' @param model Model to calculate R_tau
#' @param n Number of individuals at each observation
#' @param data Matrix of observation
#' @param t Observation time
#' @param ... Additional parameters passed to model
#' @param lower Lower bound for estimation
#' @param upper Upper bound for estimation
#' @param optimizer Optimizer to use. Default NULL(which means @seealso [optim()]
#'   will be used). If a user-defined optimizer is applied, it is highly recommended
#'   to wrap this optimizer so that it has the same parameter format as optim().
#' @param opt_arg Control list for optimizer. Default empty list
#' @param verbose Return detailed output. Default FALSE.
#'
#' @return If verbose is TRUE return the full result from optimizer otherwise return MCLE only.
#' @export
#'
#' @examples
move.MCLE <-
  function(theta,
           model,
           data,
           t,
           ...,
           lower = 0.0,
           upper = 1.0,
           optimizer = NULL,
           opt_arg = list(),
           verbose = FALSE) {
    lent <- length(t)
    if (lent != ncol(data))
      stop("Number of observation not match(data & t)")
    n <- colSums(data != 0)
    if (lent <= 20000) {
      obs <- move.pairwise(data = data, t = t, tau = TRUE)
      return(
        move.MCLE.pair(
          theta,
          model,
          n,
          obs$m,
          obs$tau,
          ...,
          lower = lower,
          upper = upper,
          optimizer = optimizer,
          opt_arg = opt_arg,
          verbose = verbose
        )
      )
    } else {
      # TODO: undecided
    }
  }

move.CI(theta, model, data, t, ..., B = 500, cl = NULL, ncore = -1) {
  if (!is.integer(B) || B <= 0) stop("B must be a positive integer")
  if (is.null(cl)) {
    if (ncore == -1) ncore <- detectCores()
    else if (!is.integer(ncore) || ncore <= 0)
      stop("ncore must be a positive integer")
    cl <- makeCluster(ncore)
  } else if (!is(cl, 'cluster')) {
    stop('cl must be a cluster')
  }
  theta_MCLE <- parSapply(cl, 0:B, function(id) {
    if (id == 0) {
      return(move.MCLE(theta, model, data, t, ..., verbose = FALSE))
    } else {
      data_bootstrap <- move.bootstrap()
    }
  })

}

#' LIR model A. The population is constant and no migrations occur.
#'
#' @param theta here theta is the estimated LIR(R_tau) directly
#' @param tau List of time interval of lagged identification
#'
#' @return rep(theta, length(tau))
#' @export
#'
#' @examples
move.model.A <- function(theta, tau) {
  return(rep(theta, length(tau)))
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
#' @param theta (\eqn{\gamma, \beta, \alpha})
#' @param tau List of time interval of lagged identification
#'
#' @return
#' @export
#'
#' @examples
move.model.C <- function(theta, tau) {
  return(theta[1]*exp(-theta[2]*tau) + theta[3])
}
