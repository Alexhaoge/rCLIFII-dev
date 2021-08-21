# Composite Likelihood Function -------------------------------------------

#' @title Composite likelihood for LIR
#'
#' @description
#' This CL function requires a sample format of pairwise lagged identification,
#' which is a list-like data containing lagged identification of of all pairs of
#' observation. In gerneral cases we should use `LIR.CL` which is implemented
#' by C++ with time complexity of \eqn{\Theta(length(tp)^2)}. `LIR.CL.pair`
#' is recommended when total observation time is very small and pairwise lagged
#' identification is easy to pre-calculate through `LIR.pairwise()`.
#' The notation of the estimation is \eqn{\hat{R_{\tau}}}.
#'
#' @param theta Parameter to calculate \eqn{\hat{R_{\tau}}}
#' @param model Model to calculate \eqn{\hat{R_{\tau}}}. Should be a function
#'   that takes `theta`, `tau_i`(and other arguments if necessary) and return a
#'   numeric number
#' @param data Observation matrix
#' @param tp List-like observation time(1d vector)
#' @param ni Number of individuals of the former observation at each lagged time pair.
#' @param nj Number of individuals of the latter observation at each lagged time pair.
#' @param m vector of total lagged identification number for all lagged time pair
#' @param tau vector of time interval of lagged identification
#' @param mtau The maximum allowable lag time. If a lagged pair has time \eqn{\tau}
#'   greater than `mtau`, it will not be used to calculate composite likelihood.
#' @param model_args extra arguments to be passed to model to calculate \eqn{\hat{R_{\tau}}}
#'
#' @return minus Likelihood score(numeric). Inf if theta out of [0,1]
#' @export
#' @rdname CL
#'
#' @examples
LIR.CL.pair <- function(theta, model, ni, nj, m, tau, model_args = list()) {
  R_tau <- sapply(tau, function(t){model(theta=theta, tau=t, model_args)})
  nR_tau = R_tau * nj
  if (min(nR_tau)<=0) return(-9e12);
  if (max(nR_tau)>=1) return(-9e12);
  likelihood <- sum(m * log(nR_tau) + (ni-m) * log(1 - nR_tau))
  return(likelihood)
}


# Gradient of Composite Likelihood ----------------------------------------

#' @title Gradient of composite likelihood
#'
#' @param theta Parameter to calculate \eqn{\hat{R_{\tau}}}
#' @param model Model to calculate \eqn{\hat{R_{\tau}}}. Should be a function
#'   that takes theta, tau_i(and other arguments if necessary) and return a
#'   numeric number
#' @param grad Gradient of the model. Should be a function that takes theta,
#'   tau_i(and other arguments if necessary) and return a 1D
#'   vector with same length as theta
#' @param data Observation matrix
#' @param tp List-like observation time(1d vector)
#' @param ni of individuals of the former observation at each lagged time pair.
#' @param nj Number of individuals of the latter observation at each lagged time pair.
#' @param m vector of total lagged identification number for all lagged time pair
#' @param tau vector of time interval of lagged identification
#' @param model_args extra arguments to be passed to model to calculate \eqn{\hat{R_{\tau}}}
#' @param mtau The maximum allowable lag time. If a lagged pair has time \eqn{\tau}
#'   greater than `mtau`, it will not be used to calculate composite likelihood.
#' @return A 1D vector with same length as theta
#' @rdname CLgrad
#' @export
#'
#' @examples
LIR.CLgrad.pair <- function(theta, model, grad, ni, nj, m, tau, model_args = list()) {
  R_tau <- sapply(tau, function(t){model(theta, t, model_args)})
  R_tau_grad <- t(sapply(tau, function(t){grad(theta, t, model_args)}))
  if (length(theta) == 1)
    R_tau_grad <- t(R_tau_grad)
  return(colSums(R_tau_grad * (m / R_tau - (ni - m) * nj / (1 - nj * R_tau))))
}


# Hessian Matrix of Composite Likelihood ----------------------------------

#' @title Hessian matrix of composite likelihood
#'
#' @param theta Parameter to calculate \eqn{\hat{R_{\tau}}}
#' @param model Model to calculate \eqn{\hat{R_{\tau}}}. Should be a function
#'   that takes theta, tau_i(and other arguments if necessary) and return a
#'   numeric number
#' @param grad Gradient of the model. Should be a function that takes theta,
#'   tau_i(and other arguments if necessary) and return a 1D
#'   vector with same length as theta
#' @param hessian Hessian of the model. Should be a function that takes theta,
#'   tau_i(and other arguments if necessary) and return a square symmetric matrix
#'   with length(theta) rows and columns
#' @param data Observation matrix
#' @param tp List-like observation time(1d vector)
#' @param ni of individuals of the former observation at each lagged time pair.
#' @param nj Number of individuals of the latter observation at each lagged time pair.
#' @param m vector of total lagged identification number for all lagged time pair
#' @param tau vector of time interval of lagged identification
#' @param model_args extra arguments to be passed to model to calculate \eqn{\hat{R_{\tau}}}
#' @param mtau The maximum allowable lag time. If a lagged pair has time \eqn{\tau}
#'   greater than `mtau`, it will not be used to calculate composite likelihood.
#'
#' @return A square symmetric matrix with length(theta) rows and columns.
#' @rdname CLhessian
#' @export
#'
#' @examples
LIR.CLhessian.pair <- function(theta, model, grad, hessian, ni, nj, m, tau, model_args = list()) {
  R_tau <- sapply(tau, function(t) {model(theta, t, model_args)})
  R_tau_grad <- sapply(tau, function(t) {grad(theta, t, model_args)}, simplify = 'matrix')
  if (length(theta) == 1)
    R_tau_grad <- t(R_tau_grad)
  d1 <- dim(R_tau_grad)
  R_tau_grad2 <-
    array(apply(R_tau_grad, 2, function(t) {t %*% t(t)}), dim = c(d1[1], d1[1], d1[2])) *
    rep((ni - m) * nj^2 / (1 - nj * R_tau)^2 + m / R_tau^2, each=d1[1]*d1[1])
  R_tau_hessian <-
    sapply(tau, function(t) {hessian(theta, t, model_args)}, simplify = 'array') *
    rep((ni - m) * nj / (1 - nj * R_tau) - m / R_tau, each=d1[1]*d1[1])
  if (d1[1] == 1) R_tau_hessian <- array(R_tau_hessian, dim = c(d1[1], d1[1], d1[2]))
  return(-rowSums(R_tau_grad2 + R_tau_hessian, dims = 2))
}


# Maximal Composite Likelihood --------------------------------------------

#' @rdname MCLE
#' @export
LIR.MCLE.pair <-
  function(theta,
           model,
           ni,
           nj,
           m,
           tau,
           model_args = list(),
           lower = 0.0,
           upper = 1.0,
           optimizer = NULL,
           opt_arg = list(),
           verbose = FALSE) {
    clf <- function(theta_) { -LIR.CL.pair(theta_, model, ni, nj, m, tau, model_args) }
    if (base::is.null(optimizer)) {
      if (length(theta) == 1)
        opt_res <-
          stats::optim(par = theta, fn = clf, method = 'Brent',
            lower = lower, upper = upper, control = opt_arg)
      else
        opt_res <- stats::optim(par = theta, fn = clf, control = opt_arg)
    } else {
      opt_res <- optimizer(theta, clf, lower = lower, upper = upper, control = opt_arg)
    }
    if(verbose) return(opt_res)
    else return(opt_res$par)
  }

#' @name LIR.MCLE
#' @title Maximal Composite Likelihood Estimate
#'
#' @description
#' `LIR.MCLE` is for observation matrix input. If number of observation(length(t))
#' is less than 20000, pairwise data will be pre-calculated and `LIR.CL.pair`
#' will be used. Otherwise pairwise lagged identification will be calculated
#' separately at every iteration.
#' `LIR.MCLE.pair` is for pairwise list data input.(Deprecated)
#'
#' @param theta Initial value for parameters to estimate
#' @param model Model to calculate \eqn{\hat{R_{\tau}}}
#' @param data Observation matrix
#' @param tp List-like observation time(1d vector)
#' @param model_args Additional parameters passed to model
#' @param mtau The maximum allowable lag time. If a lagged pair has time \eqn{\tau}
#'   greater than `mtau`, it will not be used to calculate composite likelihood.
#' @param lower Lower bound for estimation
#' @param upper Upper bound for estimation
#' @param optimizer Optimizer to use. Default NULL(which means @seealso [optim()]
#' will be used). If a user-defined optimizer is applied, it is highly recommended
#' to wrap this optimizer so that it has the same parameter format as optim().
#' @param opt_arg Control list for optimizer. Default empty list
#' @param verbose Return detailed output. Default FALSE.
#'
#' @return If verbose is TRUE return the full result from optimizer otherwise return MCLE only.
#' @export
#' @rdname MCLE
#'
#' @examples
#' # Example of MCLE
#' # Set observation time
#' t <- c(1:5, 51:55, 101:105, 501:505, 601:605)
#' # Generate observation matrix with model C
#' data <- move.simulate.C(300, 100, 605, 40, t, 0.08, 0.04)
#' # MCLE
#' LIR.MCLE(rep(0.001, 3), LIR.model.C, data, t)
#' # [1] 0.007119308 0.212260232 0.003582160
#'
LIR.MCLE <-
  function(theta,
           model,
           data,
           tp,
           model_args = list(),
           mtau = Inf,
           lower = 0.0,
           upper = 1.0,
           optimizer = NULL,
           opt_arg = list(),
           verbose = FALSE) {
    if (length(tp) != ncol(data))
      stop("Number of observation not match")
    clf <- function(theta_) {-LIR.CL(theta, model, data, tp, model_args, mtau)}
    if (base::is.null(optimizer)) {
      if (length(theta) == 1)
        opt_res <- stats::optim(par = theta, fn = clf, method = 'Brent',
            lower = lower, upper = upper, control = opt_arg)
      else
        opt_res <- stats::optim(par = theta, fn = clf, control = opt_arg)
    } else {
      opt_res <- optimizer(theta, clf, lower = lower, upper = upper, control = opt_arg)
    }
    if(verbose) return(opt_res)
    else return(opt_res$par)
  }


# Confidence Interval -----------------------------------------------------

#' Confidence interval for LIR MCLE
#'
#' @description
#' Estimate of LIR using MCLE is asymptotically normal so here bootstrap is applied to estimate the variance.
#'
#' @param theta initial value of estimate for MCLE iteration
#' @param model Model to calculate \eqn{\hat{R_{\tau}}}
#' @param data Observation matrix
#' @param tp List-like observation time(1d vector)
#' @param ... Additional parameters passed to the optimizer of MCLE
#' @param model_args Additional parameters passed to model
#' @param mtau The maximum allowable lag time. If a lagged pair has time \eqn{\tau}
#'   greater than `mtau`, it will not be used to calculate composite likelihood.
#' @param B Bootstrap's repeat sampling times
#' @param cl Cluster to use, Default NULL. If NULL, a new cluster will be created by @seealso [makeCluster()]
#' @param ncores Number of processors to use. Default -1(which means all available cores).
#'   This argument will be suppressed if cl is not NULL.
#' @param alpha Confidence level. Default 0.05.
#'
#' @return If length(theta) equals 1, return a 2 element vector representing the CI.
#'   Otherwise return a matrix. The first row of the matrix is the lower bound
#'   and the second row is the upper bound.
#' @export
#'
#' @examples
#' # Example of confidence interval of MCLE
#' # Set observation time
#' tp <- c(1:5, 51:55, 101:105, 501:505, 601:605)
#' # Generate observation matrix with model C
#' data <- LIR.simulate.C(300, 100, 40, tp, 0.08, 0.04)
#' # CI
#' LIR.CI(c(0.001, 0.001, 0.001), LIR.model.C, data, tp, B = 10)
#' #              [,1]       [,2]        [,3]
#' # lower 0.005654683 0.08993638 0.003217763
#' # upper 0.008583934 0.33458409 0.003946556
#'
LIR.CI <-
  function(theta,
           model,
           data,
           tp,
           ...,
           model_args = list(),
           mtau = Inf,
           B = 500,
           cl = NULL,
           ncores = -1,
           alpha = 0.05) {
    B <- as.integer(B)
    if (B <= 1)
      stop("B must be a positive integer bigger than 1")
    clusterNotGiven <- FALSE
    if (base::is.null(cl)) {
      clusterNotGiven <- TRUE
      if (ncores == -1)
        ncores <- parallel::detectCores()
      else if (!base::is.integer(ncores) || ncores <= 0)
        stop("ncores must be a positive integer")
      cl <- parallel::makeCluster(ncores)
    } else if (!base::is(cl, 'cluster')) {
      stop('cl must be a cluster')
    }
    theta <- LIR.MCLE(theta, model, data, tp, ..., model_args=model_args, mtau=mtau, verbose = FALSE)
    theta_boot <- parallel::parSapply(cl, 1:B, function(id) {
      data_bootstrap <- LIR.bootstrap(data)
      return(LIR.MCLE(theta, model, data_bootstrap, tp, ..., model_args=model_args, mtau=mtau, verbose = FALSE))
    })
    if (clusterNotGiven) parallel::stopCluster(cl)
    if (!base::is.numeric(alpha) || alpha <= 0 || alpha >= 1)
      stop("Confidence level should between 0 and 1(exclusive), usually a small value like 0.05")
    z <- stats::qnorm(1 - alpha / 2)
    if (length(theta) == 1) {
      std <- stats::var(theta_boot)
      return(c(theta - z * std, theta + z * std))
    } else {
      std <- sqrt(diag(stats::var(t(theta_boot))))
      ci <- rbind(theta - z * std, theta + z * std)
      rownames(ci) <- c('lower', 'upper')
      return(ci)
    }
  }

