LIR.XIC.pair <- function (theta, model, ni, nj, m, tau, ..., MCLE = FALSE) {
  if (MCLE) {
    return(-2 * LIR.CL.pair(theta, model, ni, nj, m, tau, ...))
  } else {
    mcle <- LIR.MCLE.pair(theta, model, ni, nj, m, tau, ..., verbose=TRUE)
    return(-2 * mcle$value)
  }
}

LIR.AIC.pair <- function(theta, model, ni, nj, m, tau, ..., MCLE = FALSE) {
  return(LIR.XIC.pair(theta, model, ni, nj, m, tau, ..., MCLE = MCLE) + 2 * length(theta))
}

LIR.BIC.pair <- function(theta, model, ni, nj, m, tau, `T`, ..., MCLE = FALSE) {
  return(LIR.XIC.pair(theta, model, ni, nj, m, tau, ..., MCLE = MCLE) + length(theta) * log(`T`))
}

LIR.QAIC.pair <- function(theta, model, ni, nj, m, tau, ..., MCLE = FALSE) {
  cl <- LIR.XIC.pair(theta, model, ni, nj, m, tau, ..., MCLE = MCLE)
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

#' CLIC
#'
#' @param theta
#' @param model
#' @param grad
#' @param hessian
#' @param data
#' @param t
#' @param ...
#' @param MCLE
#' @param B
#' @param cl
#' @param ncore
#' @param `T` Total study time. Required for CLICb
#'
#' @rdname CLIC
LIR.CLICbase <-
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
      cl <-
        -2 * LIR.CL.pair(theta, model, lbs$ni, obs$nj, obs$m, obs$tau, ...)
    } else {
      mcle <-
        LIR.MCLE.pair(theta, model, lbs$ni, obs$nj, obs$m, obs$tau, ..., verbose=TRUE)
      cl <- -2 * mcle$value
      theta <- mcle$par
    }
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
    theta_MCLE <- parSapply(cl, 1:B, function(id) {
      data_bootstrap <- LIR.bootstrap(data)
      return(LIR.MCLE(theta, model, data_bootstrap, t, ..., verbose = FALSE))
    })
    if (clusterNotGiven)
      stopCluster(cl)
    H <-
      LIR.CLhessian.pair(theta, model, grad, hessian, lbs$ni, obs$nj, obs$m, obs$tau, ...)
    varBoot <- var(t(theta_MCLE))
    return(c(cl, sum(diag(H %*% varBoot))))
  }

#' CLICa
#' @rdname CLIC
LIR.CLICa <-
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
    tmp <-
      LIR.CLICbase(
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

#' CLICb
#' @rdname CLIC
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
      LIR.CLICbase(
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

LIR.modelSelect <-
  function(theta,
           data,
           t,
           model_list = c('A', 'B', 'C'),
           grad_list = c('A', 'B', 'C'),
           hessian_list = c('A', 'B', 'C'),
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
    fun_list <- list()
    for (i in 1:lenf) {
      if (model_list[i] == 'A')
        fun_list[i] <- LIR.model.A
      else if (model_list[i] == 'B')
        fun_list[i] <- LIR.model.B
      else if (model_list[i] == 'C')
        fun_list[i] <- LIR.model.C
      else if (is.function(model_list[i]))
        fun_list[i] <- model_list[i]
      else
        stop('Only "A", "B", "C" or function is valid to be in `model_list`')
    }
    score <- c()
    if (criterion %in% c('AIC', 'BIC', 'QAIC')) {
      obs <- LIR.pairwise(data, t, tau = TRUE, ret_n = TRUE)
      if (criterion == 'AIC') {
        score <-
          sapply(fun_list, function(f) {
            LIR.AIC.pair(theta, f, obs$ni, obs$nj, obs$m, obs$tau, ..., MCLE)
          })
      } else if (criterion == 'BIC') {
        score <-
          sapply(fun_list, function(f) {
            LIR.BIC.pair(theta, f, obs$ni, obs$nj, obs$m, obs$tau, `T`, ..., MCLE)
          })
      } else {
          score <-
            sapply(fun_list, function(f) {
              LIR.QAIC.pair(theta, f, obs$ni, obs$nj, obs$m, obs$tau, ..., MCLE)
            })
        }
    } else if (criterion %in% c('CLICa', 'CLICb')) {
      if (lenf != length(grad_list) || lenf != length(hessian_list))
        stop("Gradient list and Hessian list must be equal length with model list")
      for (i in 1:lenf) {
        if (grad_list[i] == 'A')
          grad_list[i] <- LIR.grad.A
        else if (grad_list[i] == 'B')
          grad_list[i] <- LIR.grad.B
        else if (grad_list[i] == 'C')
          grad_list[i] <- LIR.grad.C
        else if (!is.function(grad_list[i]))
          stop('Only "A", "B", "C" or function is valid to be in `grad_list`')
        if (hessian_list[i] == 'A')
          hessian_list[i] <- LIR.hessian.A
        else if (hessian_list[i] == 'B')
          hessian_list[i] <- LIR.hessian.B
        else if (hessian_list[i] == 'C')
          hessian_list[i] <- LIR.hessian.C
        else if (!is.function(hessian_list[i]))
          stop('Only "A", "B", "C" or function is valid to be in `hessian_list`')
      }
      if (criterion == 'CLICa') {
        for (i in 1:lenf) {
          score[i] <-
            LIR.CLICa(
              theta,
              fun_list[i],
              grad_list[i],
              hessian_list[i],
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
              theta,
              fun_list[i],
              grad_list[i],
              hessian_list[i],
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
    best_idx <- which.max(score)
    if (verbose) {
      if (!MCLE)
        theta <-
          LIR.MCLE(theta, fun_list[best_idx], data, t, ..., verbose = FALSE)
      return(
        list(
          best = model_list[best_idx],
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
