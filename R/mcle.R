
# Composite Likelihood Function -------------------------------------------


#' Composite likelihood for estimating LIR(Lagged Identification Rate)
#'
#' This CL function requires a sample format of pairwise lagged identification,
#' which is a list-like data containing lagged identification of of all pairs of
#' observation.
#'
#' This function is recommended when total observation time is very small and
#' pairwise lagged identification is easy to pre-calculate through
#' @seealso [LIR.pairwise()]. The notation of the estimation is \eqn{\hat{R_{tau}}}
#'
#' @param theta Parameter to calculate \eqn{\hat{R_{tau}}}
#' @param model Model to calculate \eqn{\hat{R_{tau}}}
#' @param ni Number of individuals at each observation
#' @param nj Number of individuals at each lagged observation
#' @param m vector of lagged identification for all observation pair
#' @param tau vector of time interval of lagged identification
#' @param ... extra arguments to be passed to model to calculate \eqn{\hat{R_{tau}}}
#'
#' @return Likelihood score(numeric). Inf if theta out of [0,1]
#' @export
#' @rdname CL
#' @examples
LIR.CL.pair <- function(theta, model, ni, nj, m, tau, ...) {
  R_tau <- sapply(tau, function(t){model(theta, t, ...)})
  if (min(R_tau * nj)<=0) return(Inf);
  if (max(R_tau * nj)>=1) return(Inf);
  likelihood <- sum(m * log(R_tau*nj) + (ni-m) * log(1 - R_tau*nj))
  return(-likelihood)
}

#' #' @rdname CL
#' LIR.CL <- function(theta, model, n, data, t, ...) {
#'
#' }


# Gradient of Composite Likelihood ----------------------------------------

#' Title
#'
#' @param theta
#' @param model
#' @param grad
#' @param ni
#' @param nj
#' @param m
#' @param tau
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
LIR.CLgrad.pair <- function(theta, model, grad, ni, nj, m, tau, ...) {
  R_tau <- sapply(tau, function(t){model(theta, t, ...)})
  R_tau_grad <- t(sapply(tau, function(t){grad(theta, t, ...)}))
  return(-colSums(R_tau_grad * ((ni - m) * nj / (1 - nj * R_tau) - mj / R_tau)))
}


# Hessian Matrix of Composite Likelihood ----------------------------------

#' Title
#'
#' @param theta
#' @param model
#' @param grad
#' @param hessian
#' @param ni
#' @param nj
#' @param m
#' @param tau
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
LIR.CLhessian.pair <- function(theta, model, grad, hessian, ni, nj, m, tau, ...) {
  R_tau <- sapply(tau, function(t) {model(theta, t, ...)})
  R_tau_grad <- sapply(tau, function(t) {grad(theta, t, ...)})
  d1 <- dim(R_tau_grad)
  R_tau_grad2 <-
    array(apply(ss, 2, function(t) {t %*% t(t)}), dim = c(d1[1], d1[1], d1[2])) *
    rep((ni - m) * nj^2 / (1 - nj * R_tau)^2 + m / R_tau^2, each=d1[1]*d[1])
  R_tau_hessian <-
    sapply(tau, function(t) {hessian(theta, t, ...)}, simplify = 'array') *
    rep((ni - m) * nj / (1 - nj * R_tau) - m / R_tau, each=d1[1]*d[1])
  return(-rowSums(R_tau_grad2 + R_tau_hessian, dims = 2))
}


# Maximal Composite Likelihood --------------------------------------------

#' Maximal Composite Likelihood Estimate for pairwise data
#'
#' @param theta Initial value for parameters to estimate
#' @param model Model to calculate \eqn{\hat{R_{tau}}}
#' @param ni Number of individuals at each observation
#' @param nj Number of individuals at each lagged observation
#' @param m Vector of lagged identification for all observation pair
#' @param tau Vector of time interval of lagged identification
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
LIR.MCLE.pair <-
  function(theta,
           model,
           ni,
           nj,
           m,
           tau,
           ...,
           lower = 0.0,
           upper = 1.0,
           optimizer = NULL,
           opt_arg = list(),
           verbose = FALSE) {
    clf <-
      function(theta_) {
        LIR.CL.pair(
          theta_,
          model = model,
          ni = ni,
          nj = nj,
          m = m,
          tau = tau,
          ...
        )
      }
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
      opt_res <-
        optimizer(theta,
                  clf,
                  lower = lower,
                  upper = upper,
                  control = opt_arg)
    }
    return(ifelse(verbose, opt_res, opt_res$par))
  }

#' Maximal Composite Likelihood Estimate
#' If number of observation(length(t)) is less than 20000, pairwise MCLE
#' will be used. Otherwise pairwise lagged identification will be
#' calculated seperately at every iteration.
#'
#' @param theta Initial value for parameters to estimate
#' @param model Model to calculate \eqn{\hat{R_{tau}}}
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
LIR.MCLE <-
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
      obs <- LIR.pairwise(data = data, t = t, tau = TRUE)
      return(
        LIR.MCLE.pair(
          theta,
          model,
          obs$ni,
          obs$nj,
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
      stop('Too much observation')
      # TODO: undecided
    }
  }


# Confidence Interval -----------------------------------------------------

#' Confidence interval for LIR MCLE
#' Estimate of LIR using MCLE is asymptotically normal so here bootstrap is applied to estimate the variance.
#'
#' @param theta initial value of estimate for MCLE iteration
#' @param model Model to calculate \eqn{\hat{R_{tau}}}
#' @param data Matrix of observation
#' @param t Observation time
#' @param ...
#' @param B Bootstrap's repeat sampling times
#' @param cl Cluster to use, Default NULL. If NULL, a new cluster will be created by @seealso [makeCluster()]
#' @param ncore Number of processors to use. Default -1(which means all available cores).
#'   This argument will be suppressed if cl is not NULL.
#' @param alpha Confidence level. Default 0.05.
#'
#' @return
#' @export
#'
#' @examples
LIR.CI <-
  function(theta,
           model,
           data,
           t,
           ...,
           B = 500,
           cl = NULL,
           ncore = -1,
           alpha = 0.05) {
    if (!is.integer(B) || B <= 0)
      stop("B must be a positive integer")
    clusterNotGiven <- FALSE
    if (is.null(cl)) {
      clusterNotGiven <- TRUE
      if (ncore == -1)
        ncore <- detectCores()
      else if (!is.integer(ncore) || ncore <= 0)
        stop("ncore must be a positive integer")
      cl <- makeCluster(ncore)
    } else if (!is(cl, 'cluster')) {
      stop('cl must be a cluster')
    }
    theta_MCLE <- parSapply(cl, 0:B, function(id) {
      if (id == 0) {
        return(LIR.MCLE(theta, model, data, t, ..., verbose = FALSE))
      } else {
        data_bootstrap <- LIR.bootstrap(data)
        return(LIR.MCLE(theta, model, data_bootstrap, t, ..., verbose = FALSE))
      }
    })
    if (clusterNotGiven) stopCluster(cl)
    if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1)
      stop("Confidence level should between 0 and 1(exclusive), usually a small value like 0.05")
    z <- qnorm(1 - alpha / 2)
    if (length(theta) == 1) {
      std <- var(theta_MCLE[-1])
      return(c(theta_MCLE[1] - z * std, theta_MCLE[1] + z * std))
    } else {
      std <- sqrt(diag(var(t(theta_MCLE[,-1]))))
      ci <- rbind(theta_MCLE[, 1] - z * std, theta_MCLE[, 1] + z * std)
      colnames(ci) <- c('lower', 'upper')
      return(ci)
    }
  }

