
#' @name LIR.simulate.A/B/C
#' @title Simulate animal movement identification data with model A/B/C
#'
#' @param Z Total population of the whole area.(Only for model C)
#' @param N Population within the study area.(For model C, the initial population)
#' @param n Number of identification in each observation
#' @param tp Time of each observation if tp is a list or array, otherwise the number
#'   of observation if t is an integer. Model B and C require tp to be a list.
#' @param lambda move-out rate
#' @param mu move-in rate
#' @param seed random seed
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
#' @return matrix of length(t) columns. The number of rows is N for model A/B,
#' Z for model C.
#' @export
#' @rdname simulate
#'
#' @examples
#' # Example of data simulation
#' # Set observation time
#' t <- c(1:5, 51:55, 101:105, 501:505, 601:605)
#' # Generate observation matrix with model C
#' data <- move.simulate.C(T=300, N=100, n=40, tp=t, lambda=0.08, mu=0.04)
#' dim(data)
#' # [1] 300  25
#'
LIR.simulate.A <- function(N, n, tp, seed = NULL) {
  set.seed(seed)
  lent <- ifelse(length(tp) == 1, tp, length(tp))
  if (length(n) == 1)
    data <- replicate(lent, (sampling::srswor(n, N) > 0) * 1)
  else if (length(n) == lent)
    data <- sapply(n, FUN=function(ni){(sampling::srswor(ni, N) > 0) * 1})
  else stop("Dimension of n does not equal t")
  return(data)
}

#' @rdname simulate
#' @export
LIR.simulate.B <- function(N, n, tp, lambda, seed = NULL) {
  set.seed(seed)
  lent <- length(tp)
  if (lent == 1)
    stop("For model B, t must be a list")
  else if (length(n) > 1 && length(n) != lent)
    stop("Dimension of n does not equal t")
  if (length(n == 1)) n = rep(n, lent)
  data <- matrix(nrow = N, ncol = lent)
  pop <- rep(1, N)
  j <- 1
  tp <- sort(tp)
  if (1 == tp[j]) {
    data[, 1] <- pop * (sampling::srswor(n[1], N) > 0)
    j <- j + 1
  }
  tT <- tp[lent]
  for (i in 2:tT) {
    unif <- stats::runif(N, 0, 1)
    pop <- (unif <= lambda) + pop
    if (tp[j] == i) {
      data[, j] <- pop * (sampling::srswor(n[j], N) > 0)
      j <- j + 1
      if (j > lent) break
    }
  }
  return(data)
}

#' @rdname simulate
#' @export
LIR.simulate.C <- function(Z, N, n, tp, lambda, mu, seed = NULL) {
  set.seed(seed)
  lent <- length(tp)
  if (lent == 1)
    stop("For model C, t must be a list")
  if (length(n) > 1 && length(n) != lent)
    stop("Dimension of n does not equal t")
  if (length(n == 1)) n = rep(n, lent)
  if (lambda / mu != Z/N - 1)
    warning("Lambda, mu given cannot gurantee that expectancy of
            population in invesgated area remains N")
  data <- matrix(0, nrow = Z, ncol = lent)
  pop_in <- seq(N)
  pop_out <- seq(N+1, Z)
  j <- 1
  tp <- sort(tp)
  if (1 == tp[j]) {
    data[sample(pop_in, n[1]), 1] <- 1
    j <- j + 1
  }
  tT <- tp[lent]
  for (i in 2:tT) {
    unif_in <- stats::runif(length(pop_in))
    unif_out <- stats::runif(length(pop_out))
    pop_in_new <- c(
      pop_in[which(unif_in > lambda)],
      pop_out[which(unif_out <= mu)]
    )
    pop_out_new <- c(
      pop_in[which(unif_in <= lambda)],
      pop_out[which(unif_out > mu)]
    )
    pop_in <- pop_in_new
    pop_out <- pop_out_new
    if (i == tp[j]) {
      data[sample(pop_in, n[j]), j] <- 1
      j <- j + 1
      if (j > lent) break
    }
  }
  return(data)
}

#' @title Calculate lagged identification of each observation pair.
#' @description
#' Calculate lagged identification of each observation pair. Number of observation is guessed by the column number of input matrix
#'
#' @param data Matrix of identification with row number equal to population
#'   and column number equal to observation.
#' @param tp List like time of observation. Default NULL, must given if tau is true
#' @param require_tau If true, return tau
#' @param require_n If true, return n
#'
#' @return List of lagged identification number for each observation pair.
#' @export
#'
#' @examples
LIR.pairwise <- function(data, tp = NULL, require_tau = TRUE, require_n = TRUE) {
  lent <- length(tp)
  if (ncol(data) != lent)
    stop("Number of data columns does not match with length of t")
  obs  <- list()
  obs[["m"]] <- unlist(
    sapply(1:(lent-1), FUN=function(i) {
      sapply((i+1):lent, FUN = function(j){
        sum((data[, i] == data[, j]) * (data[, i] != 0))
      })
    })
  )
  if (require_tau) {
    if (is.null(tp))
      stop("No observation time provided to calculate tau")
    obs[["tau"]] <- unlist(
      sapply(1:(lent-1), FUN=function(i) {
        sapply((i+1):lent, FUN = function(j){return(tp[j] - tp[i])})
      })
    )
  }
  if (require_n) {
    n <- colSums(data)
    obs[["ni"]] <- unlist(
      sapply(1:(lent-1), FUN=function(i) {
        sapply((i+1):lent, FUN = function(j){n[i]})
      })
    )
    obs[["nj"]] <- unlist(
      sapply(1:(lent-1), FUN=function(i) {
        sapply((i+1):lent, FUN = function(j){n[j]})
      })
    )
  }
  return(obs)
}


#' Non-parametric estimation of LIR
#'
#' @param data
#' @param n
#' @param tp
#'
#' @return
#' @export
#'
#' @examples
non.lir <- function(data, n, tp) {


  lent <- length(tp)
  tp <- sort(tp)
  tT <- tp[lent]

  R.m <- rep(0, tT)
  R.n <- rep(0, tT)
  mij <- c()
  nij <- c()
  tauij <- c()

  k <- 1
  for (i in 1:(lent-1)){
    for (j in (i+1):lent){
      mij[k] <- sum((data[, i] == data[, j]) * (data[, i] != 0))
      nij[k] <- n[i]*n[j]
      tauij[k] <- tp[j] - tp[i]
      R.m[tauij[k]] <- R.m[tauij[k]] + mij[k]
      R.n[tauij[k]] <- n[i]*n[j] + R.n[tauij[k]]
      k <- k + 1
    }

  }

  R.tauij <- R.m[tauij]/R.n[tauij]
  tau <- unique(tauij)
  R.tau <- R.m[tau]/R.n[tau]

  R.dat <- list(R.tau=R.tau, tau=tau, R.m=R.m, R.n=R.n,
                mij=mij, nij=nij, tauij=tauij)

  return(R.dat)
}

#' @title Bootstrap on individuals
#'
#' @description
#' Perform a single bootstrap on the given observation matrix.
#' SRSWR will be performed on individuals (row of the matrix)
#'
#' @param data Observation matrix. Each row represent an individual
#'  and each column represent an observation
#' @param seed random seed
#'
#' @return An new observation matrix. Same dims as the input data/
#' @export
#'
#' @examples
LIR.bootstrap <- function(data, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(data)
  return(rbind(data[sample(1:n, n, replace=TRUE),]))
}

# LIR.plotLIR <- function(data, Time, n, t, fun.R.tau, ..., title = NULL) {
#   # TODO: codes here is slow and need cleaning
#   R.m <- rep(0,Time)
#   R.n <- rep(0,Time)
#   m <- c()
#   tau <- c()
#   k <- 1
#   for (i in 1:(length(t)-1)){
#     for (j in (i+1):length(t)){
#       m[k] <- sum((data[, i] == data[, j]) * (data[, i] != 0))
#       tau[k] <- t[j]-t[i]
#       R.m[tau[k]] <- R.m[tau[k]] + m[k]
#       R.n[tau[k]] <- n[i]*n[j] + R.n[tau[k]]
#       k <- k + 1
#     }
#   }
#   R.tau <- R.m[tau]/R.n[tau]
#   plot(tau, R.tau,ylim=c(0,0.02) ,ylab="Lagged Identification Rates", xlab="Time lag")
#   y = sapply(seq(Time), function(x){fun.R.tau(x, ...)})
#   lines(y, col='red', lwd=1.5)
#   legend('topright', title)
# }

