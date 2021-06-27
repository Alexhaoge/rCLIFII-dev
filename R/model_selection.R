
# AIC/BIC/QAIC ------------------------------------------------------------

# Base function for AIC/BIC/QAIC. Should not export
.LIR.XIC.pair <- function (theta, model, ni, nj, m, tau, ..., MCLE = FALSE) {
  if (MCLE) {
    return(2 * LIR.CL.pair(theta, model, ni, nj, m, tau, ...))
  } else {
    mcle <- LIR.MCLE.pair(theta, model, ni, nj, m, tau, ..., verbose=TRUE)
    return(2 * mcle$value)
  }
}

#' @name  LIR.XIC.pair
#' @title Other model selection criterion(AIC/BIC/QAIC)
#'
#' Not recommended for composite likelihood model.
#'
#' @description
#' `LIR.AIC.pair` is for AIC with pairwise data.
#' `LIR.BIC.pair` is for BIC with pairwise data.
#' `LIR.QAIC.pair` is for QIC with pairwise data.
#' \deqn{AIC(\hat{\theta}_{LIR})=-2\sum_{t_i \in \tau_0}\sum_{\tau \in M}\{cl_{ij}(\hat{\theta}_{LIR})\}+2k}
#' \deqn{BIC(\hat{\theta}_{LIR})=-2\sum_{t_i \in \tau_0}\sum_{\tau \in M}\{cl_{ij}(\hat{\theta}_{LIR})\}+k ln(T)}
#' \deqn{QAIC(\hat{\theta}_{LIR})=-2\sum_{t_i \in \tau_0}\sum_{\tau \in M}\{cl_{ij}(\hat{\theta}_{LIR})\}/\hat{c}+2k}
#' Here \eqn{k} is `length(theta)`, T is total study time, \eqn{\hat{c}} is the
#' variance inflation factor estimated from the ratio of the goodness-of-fit
#' \eqn{\chi^2}-statistic to its degrees of freedom.
#'
#' @param theta Parameter to calculate \eqn{\hat{R_{\tau}}}. If param `MCLE` is TRUE,
#' this is the MCLE, otherwise the initial value for optimizer. Theta is required
#' because it is the only parameter that tells the dimension of estimating variable.
#' @param model Model to calculate \eqn{\hat{R_{\tau}}}. Should be a function
#' that takes theta, tau_i(and other arguments if necessary) and return a
#' numeric number
#' @param ni Number of individuals at each observation
#' @param nj Number of individuals at each lagged observation
#' @param m Vector of lagged identification for all observation pair
#' @param tau Vector of time interval of lagged identification
#' @param T Total study time. Required for BIC.
#' @param ... Additional parameters passed to model, or upper/lower bound and optimize
#' argument when calculating MCLE.
#' @param MCLE If TRUE param `theta` is the MCLE, otherwise `theta` will be used as
#' the initial value for optimizer. Boolean, default FALSE.
#'
#' @return score
#' @export
#' @rdname AICBICQAIC
#'
#' @examples
LIR.AIC.pair <- function(theta, model, ni, nj, m, tau, ..., MCLE = FALSE) {
  return(.LIR.XIC.pair(theta, model, ni, nj, m, tau, ..., MCLE = MCLE) + 2 * length(theta))
}

#' @rdname AICBICQAIC
#' @export
LIR.BIC.pair <- function(theta, model, ni, nj, m, tau, `T`, ..., MCLE = FALSE) {
  return(.LIR.XIC.pair(theta, model, ni, nj, m, tau, ..., MCLE = MCLE) + length(theta) * log(`T`))
}

#' @rdname AICBICQAIC
#' @export
LIR.QAIC.pair <- function(theta, model, ni, nj, m, tau, ..., MCLE = FALSE) {
  cl <- .LIR.XIC.pair(theta, model, ni, nj, m, tau, ..., MCLE = MCLE)
  tab <- table(tau)
  ltab <- tab
  df <- ltab - length(theta) - 1
  x2 <- 0
  tmp <- c()
  cnt_category <- 0
  for (i in 1:ltab) {
    cnt_category <- cnt_category + tab[i]
    append(tmp, as.numeric(names(tab[i])))
    if (cnt_category >= 6 || i == ltab) {
      m_itau <- m[tau %in% tmp]
      mean_m <- mean(m_itau)
      x2 <- x2 + sum((m_itau - mean_m)^2) / mean_m
      cnt_category <- 0
      tmp <- c()
    }
  }
  c <- x2 / df
  return(cl / c + 2 * length(theta))
}


# CLIC --------------------------------------------------------------------

# Base function for LIR.CLICa and LIR.CLICb. Should not export
.LIR.CLICbase <-
  function(theta,
           model,
           grad,
           hessian,
           data,
           t,
           ...,
           MCLE = FALSE,
           B = 500,
           cl = NULL,
           ncore = -1) {
    obs <- LIR.pairwise(data, t, tau = TRUE, n = TRUE)
    if (MCLE) {
      clh <-
        2 * LIR.CL.pair(theta, model, obs$ni, obs$nj, obs$m, obs$tau, ...)
    } else {
      mcle <-
        LIR.MCLE.pair(theta, model, obs$ni, obs$nj, obs$m, obs$tau, ..., verbose=TRUE)
      clh <- 2 * mcle$value
      theta <- mcle$par
    }
    B <- as.integer(B)
    if (B <= 1)
      stop("B must be a positive integer bigger than 1")
    clusterNotGiven <- FALSE
    if (is.null(cl)) {
      clusterNotGiven <- TRUE
      if (ncore == -1)
        ncore <- parallel::detectCores()
      else if (!is.integer(ncore) || ncore <= 0)
        stop("ncore must be a positive integer")
      cl <- parallel::makeCluster(ncore)
    } else if (!is(cl, 'cluster')) {
      stop('cl must be a cluster')
    }
    theta_MCLE <- parallel::parSapply(cl, 1:B, function(id) {
      data_bootstrap <- LIR.bootstrap(data)
      return(LIR.MCLE(theta, model, data_bootstrap, t, ..., verbose = FALSE))
    })
    if (clusterNotGiven)
      parallel::stopCluster(cl)
    H <-
      LIR.CLhessian.pair(theta, model, grad, hessian, obs$ni, obs$nj, obs$m, obs$tau, ...)
    if (length(theta) > 1)
      theta_MCLE <- t(theta_MCLE)
    varBoot <- stats::var(theta_MCLE)
    return(c(clh, sum(diag(H %*% varBoot))))
  }


#' @name LIR.CLICa/b
#' @title Composite likelihood Information Criterion
#'
#' @description
#' `LIR.CLICa` is CLICa criterion and `LIR.CLICb` is CLICb criterion.
#' \deqn{CLIC_a(\hat{\theta}_{LIR})=-2\sum_{t_i \in \tau_0}\sum_{\tau \in M}\{cl_{ij}(\hat{\theta}_{LIR})\}+2tr\{J_T(\theta_{LIR})H_T(\theta_{LIR})^{-1}\}}
#' \deqn{CLIC_b(\hat{\theta}_{LIR})=-2\sum_{t_i \in \tau_0}\sum_{\tau \in M}\{cl_{ij}(\hat{\theta}_{LIR})\}+ln(T)tr\{J_T(\theta_{LIR})H_T(\theta_{LIR})^{-1}\}}
#' \eqn{tr\{J_T(\theta_{LIR})H_T(\theta_{LIR})^{-1}\}} is difficult to calculate
#' thus replaced by \eqn{tr\{H_T(\theta_{LIR})var(\hat{\theta}_{LIR})\}}, which is more
#' robust. Bootstrap of individual is applied to estimate covariance matrix. Re-sample
#' procedure is paralleled.
#'
#' @param theta Parameter to calculate \eqn{\hat{R_{\tau}}}. If param `MCLE` is TRUE,
#' this is the MCLE, otherwise the initial value for optimizer. Theta is required
#' because it is the only parameter that tells the dimension of estimating variable.
#' @param model Model to calculate \eqn{\hat{R_{\tau}}}. Should be a function
#' that takes theta, tau_i(and other arguments if necessary) and return a
#' numeric number
#' @param grad Gradient of the model. Should be a function that takes theta,
#' tau_i(and other arguments if necessary) and return a 1D
#' vector with same length as theta
#' @param hessian Hessian of the model. Should be a function that takes theta,
#' tau_i(and other arguments if necessary) and return a square symmetric matrix
#' with length(theta) rows and columns
#' @param data Observation matrix
#' @param t List-like observation time(1d vector)
#' @param ... Additional parameters passed to model, or upper/lower bound and optimize
#' argument when calculating MCLE.
#' @param MCLE If TRUE param `theta` is the MCLE, otherwise `theta` will be used as
#' the initial value for optimizer. Boolean, default FALSE.
#' @param B Bootstrap's repeat sampling times
#' @param cl Cluster to use, Default NULL. If NULL, a new cluster will be created by @seealso [makeCluster()]
#' @param ncore Number of processors to use. Default -1(which means all available cores).
#' This argument will be suppressed if cl is not NULL.
#' @param `T` Total study time. Required for CLICb.
#'
#' @return CLIC score
#' @export
#' @rdname CLIC
#'
#' @example
LIR.CLICa <-
  function(theta, model, grad, hessian, data, t, ..., MCLE = FALSE, B = 500, cl = NULL, ncore = -1) {
    tmp <-
      .LIR.CLICbase(
        theta,
        model,
        grad,
        hessian,
        data,
        t,
        ...,
        MCLE = MCLE,
        B = B,
        cl = cl,
        ncore = ncore
      )
    return(tmp[1] + 2 * tmp[2])
  }

#' @rdname CLIC
#' @export
LIR.CLICb <-
  function(theta,
           model,
           grad,
           hessian,
           data,
           t,
           `T`,
           ...,
           MCLE = FALSE,
           B = 500,
           cl = NULL,
           ncore = -1) {
    tmp <-
      .LIR.CLICbase(
        theta,
        model,
        grad,
        hessian,
        data,
        t,
        ...,
        MCLE = MCLE,
        B = B,
        cl = cl,
        ncore = ncore
      )
    return(tmp[1] + log(`T`) * tmp[2])
  }


# LIR Model Selection -----------------------------------------------------

#' Model Selection for lagged identification rate model.
#'
#' @description
#' List of model should be provided and a specific criterion among
#' AIC/BIC/QAIC/CLICa/CLICb should be assigned. This function will calculate
#' score of each model with the given criterion and select the model with highest
#' score as the best one.
#'
#' @param theta_list List of Parameter to calculate \eqn{\hat{R_{\tau}}}
#' correspond to `model_list`. If param `MCLE` is TRUE,
#' this is the MCLE, otherwise the initial value for optimizer. Theta is required
#' because it is the only parameter that tells the dimension of estimating variable.
#' @param data Observation matrix
#' @param t List-like observation time(1d vector)
#' @param model_list List of candidate models, default list('A', 'B', 'C').
#' Only `as.list(function)` or 'A'/'B'/'C' is accepted.'A', 'B', 'C' are
#' package's built-in models and will be replaced by `LIR.model.A`,
#' `LIR.model.B`, `LIR.model.B`, respectively.
#' @param grad_list List of gradients correspond to candidate models.
#' Only function or 'A'/'B'/'C' is accepted.'A', 'B', 'C' are package's built-in
#' models and will be replaced by `LIR.grad.A`, `LIR.grad.B`, `LIR.grad.B`, respectively.
#' Only required for CLICa/CLICb.
#' @param hessian_list List of Hessian correspond to candidate models. Default list('A', 'B', 'C').
#' Only function or 'A'/'B'/'C' is accepted.'A', 'B', 'C' will be replaced by
#' `LIR.hessian.A`, `LIR.hessian.B`, `LIR.hessian.B`, respectively. Only required for CLICa/CLICb.
#' @param criterion Model selection criterion. Should be one of 'AIC', 'BIC',
#' 'QAIC', 'CLICa' and 'CLICb'. Default list('A', 'B', 'C').
#' @param ... Additional parameters passed to model, or upper/lower bound and optimize
#' argument when calculating MCLE.
#' @param MCLE If TRUE param `theta` is the MCLE, otherwise `theta` will be used as
#' the initial value for optimizer. Boolean, default FALSE.
#' @param `T` Total study time. Required for CLICb and BIC.
#' @param B Bootstrap's repeat sampling times
#' @param cl Cluster to use, Default NULL. If NULL, a new cluster will be created by @seealso [makeCluster()]
#' @param ncore Number of processors to use. Default -1(which means all available cores).
#'   This argument will be suppressed if cl is not NULL.
#' @param verbose Whether return verbose output. Default TRUE.
#'
#' @note
#' User-defined functions in `model_list`, `grad_list`, `hessian_list` must be passed
#' as a list using as.list(function) and when iterate through model_list,
#' list format user-defined function will be convert to function by as.function(list)
#'
#' @return If `verbose` is FALSE return the score list of each model only. Otherwise
#' return a list.
#' \describe{
#'   \item{best}{Best model name}
#'   \item{best_model}{Best model function}
#'   \item{theta}{MCLE param to the best model}
#'   \item{score}{List of score to each model in `model_list`}
#'   \item{model_list}{Candidate model list}
#'   \item{criterion}{Criterion of model selection, characters}
#' }
#' @export
#'
#' @examples
#' # Example of model selection
#' # Set observation time
#' t <- c(1:5, 51:55, 101:105, 501:505, 601:605)
#' # Generate observation matrix with model C
#' data <- move.simulate.C(300, 100, 605, 40, t, 0.08, 0.04)
#' # Model Selection
#' # Initial value of theta and boundary(when length(theta)==1) are very important
#' ms <- LIR.modelSelect(data, t, B=4, MCLE=FALSE, lower=0.0005, upper=0.09)
#' ms$score
#' # [1] 1.800000e+13 1.099621e+04 1.074243e+04
#' ms$best
#' # [1] "C"
#' ms$theta
#' # [1] 0.007119308 0.212260232 0.003582160
#'
LIR.modelSelect <-
  function(data,
           t,
           theta_list = list(rep(0.1, 1), rep(0.001, 2), rep(0.001, 3)),
           model_list = list('A', 'B', 'C'),
           grad_list = list('A', 'B', 'C'),
           hessian_list = list('A', 'B', 'C'),
           criterion = 'CLICa',
           ...,
           MCLE = FALSE,
           `T` = NULL,
           B = 500,
           cl = NULL,
           ncore = -1,
           verbose = TRUE) {
    if ((criterion == 'CLICb' || criterion == 'BIC') && is.null(`T`))
      stop("Total time `T` is needed for BIC/CLICb")
    lenf <- length(model_list)
    if (lenf != length(theta_list))
      stop("Theta list should have same length with model_list")
    fun_list <- list()
    for (i in 1:lenf) {
      if (model_list[i] == 'A')
        fun_list[[i]] <- as.list(LIR.model.A)
      else if (model_list[i] == 'B')
        fun_list[[i]] <- as.list(LIR.model.B)
      else if (model_list[i] == 'C')
        fun_list[[i]] <- as.list(LIR.model.C)
      else if ('try-error' %in% class(try(as.function(grad_list[i]))))
        stop('Only "A", "B", "C" or as.list(function) is valid to be in `model_list`')
      else
        fun_list[i] <- model_list[i]
    }
    score <- c()
    if (criterion %in% c('AIC', 'BIC', 'QAIC')) {
      obs <- LIR.pairwise(data, t, tau = TRUE, n = TRUE)
      if (criterion == 'AIC') {
        for(i in 1:lenf) {
          score[i] <- LIR.AIC.pair(theta_list[[i]], as.function(fun_list[[i]]), obs$ni, obs$nj, obs$m, obs$tau, ..., MCLE)
        }
      } else if (criterion == 'BIC') {
        for(i in 1:lenf) {
          score[i] <- LIR.BIC.pair(theta_list[[i]], as.function(fun_list[[i]]), obs$ni, obs$nj, obs$m, obs$tau, `T`, ..., MCLE)
        }
      } else {
        for(i in 1:lenf) {
          score[i] <- LIR.QAIC.pair(theta_list[[i]], as.function(fun_list[[i]]), obs$ni, obs$nj, obs$m, obs$tau, ..., MCLE)
        }
      }
    } else if (criterion %in% c('CLICa', 'CLICb')) {
      if (lenf != length(grad_list) || lenf != length(hessian_list))
        stop("Gradient list and Hessian list must be equal length with model list")
      for (i in 1:lenf) {
        if (grad_list[i] == 'A')
          grad_list[[i]] <- as.list(LIR.grad.A)
        else if (grad_list[i] == 'B')
          grad_list[[i]] <- as.list(LIR.grad.B)
        else if (grad_list[i] == 'C')
          grad_list[[i]] <- as.list(LIR.grad.C)
        else if ('try-error' %in% class(try(as.function(grad_list[i]))))
          stop('Only "A", "B", "C" or function is valid to be in `grad_list`')
        if (hessian_list[i] == 'A')
          hessian_list[[i]] <- as.list(LIR.hessian.A)
        else if (hessian_list[i] == 'B')
          hessian_list[[i]] <- as.list(LIR.hessian.B)
        else if (hessian_list[i] == 'C')
          hessian_list[[i]] <- as.list(LIR.hessian.C)
        else if ('try-error' %in% class(try(as.function(hessian_list[i]))))
          stop('Only "A", "B", "C" or function is valid to be in `hessian_list`')
      }
      if (criterion == 'CLICa') {
        for (i in 1:lenf) {
          score[i] <-
            LIR.CLICa(
              theta_list[[i]],
              as.function(fun_list[[i]]),
              as.function(grad_list[[i]]),
              as.function(hessian_list[[i]]),
              data,
              t,
              ...,
              MCLE = MCLE,
              B = B,
              cl = cl,
              ncore = ncore
            )
        }
      } else {
        for (i in 1:lenf) {
          score[i] <-
            LIR.CLICb(
              theta_list[[i]],
              as.function(fun_list[[i]]),
              as.function(grad_list[[i]]),
              as.function(hessian_list[[i]]),
              data,
              t,
              `T`,
              ...,
              MCLE = MCLE,
              B = B,
              cl = cl,
              ncore = ncore
            )
        }
      }
    } else {
      stop('Unsupport criterion, only AIC,BIC, QAIC, CLICa, CLICb is allowed')
    }
    best_idx <- which.min(score)
    if (verbose) {
      if (!MCLE)
        theta <-
          LIR.MCLE(theta_list[[best_idx]], as.function(fun_list[[best_idx]]), data, t, ..., verbose = FALSE)
      return(
        list(
          best = model_list[[best_idx]],
          best_model = fun_list[[best_idx]],
          theta = theta,
          score = score,
          models = model_list,
          criterion = criterion
        )
      )
    }
    else
      return(score)
  }
