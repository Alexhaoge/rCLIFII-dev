
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
#' tp <- c(1:5, 51:55, 101:105, 501:505, 601:605)
#' # Generate observation matrix with model C
#' data <- LIR.simulate.C(Z=300, N=100, n=40, tp=tp, lambda=0.08, mu=0.04)
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
#' @param mtau The maximum allowable lag time. If a lagged pair has time \eqn{\tau}
#'   greater than `mtau`, it will not be calculated.
#'
#' @return List of lagged identification number for each observation pair.
#' @export
#'
#' @examples
LIR.pairwise <- function(data, tp = NULL, mtau = Inf) {
  lent <- length(tp)
  tp <- sort(tp)
  if (ncol(data) != lent)
    stop("Number of data columns does not match with length of t")
  m <- c()
  ni <- c()
  nj <- c()
  tau <- c()
  k <- 1
  for (i in 1:(lent-1)){
    for (j in (i+1):lent)
      if (tp[j] - tp[i] <= mtau){
        m[k] <- sum((data[, i] == data[, j]) * (data[, i] != 0))
        tau[k] <- tp[j] - tp[i]
        ni[k] <- n[i]
        nj[k] <- n[j]
        k <- k + 1
      }
  }
  return(list(m=m, ni=ni, nj=nj, tau=tau))
}


#' Non-parametric estimation of LIR
#'
#' @param data
#' @param tp
#'
#' @return
#' @export
#'
#' @examples
non.lir <- function(data, tp, mtau = Inf) {
  n <- colSums(data)
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

