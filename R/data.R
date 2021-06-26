
#' Simulate animal movement identification data with model A
#' Population is constant and no migration.
#'
#' @param N population
#' @param n number of observed individual at each moment. List or integer
#' @param t Time of each observation if t is a list or array, the number
#'   of observation if t is an integer.
#' @param seed random seed
#'
#' @return 0-1 identification matrix of N rows and length_t(or t) columns
#' @export
#'
#' @examples
LIR.simulate.A <- function(N, n, t, seed = NULL) {
  set.seed(seed)
  lent <- ifelse(length(t) == 1, t, length(t))
  if (length(n) == 1)
    data <- replicate(lent, ifelse(srswor(n, N) > 0, 1, 0))
  else if (length(n) == lent)
    data <- sapply(n, FUN=function(ni){ifelse(srswor(ni, N) > 0, 1, 0)})
  else stop("Dimension of n does not equal t")
  return(data)
}

#' Simulate animal movement identification data with model B
#'
#' @param N Population with the study area
#' @param T total time
#' @param n Number of identification in each observation
#' @param t Time of each observation if t is a list or array, the number
#'   of observation if t is an integer.
#' @param lambda migration rate
#' @param seed random seed
#'
#' @return matrix of N rows and length_t(or t) columns
#' @export
#'
#' @examples
LIR.simulate.B <- function(N, T, n, t, lambda, seed = NULL) {
  set.seed(seed)
  lent <- length(t)
  if (lent == 1) warning("For model B, t had better be a list")
  if (length(n) > 1 && length(n) != lent)
    stop("Dimension of n does not equal t")
  if (length(n == 1)) n = rep(n, lent)
  data <- matrix(nrow = N, ncol = lent)
  pop <- seq(N)
  j <- 1
  t <- sort(t)
  if (1 == t[j]) {
    data[, 1] <- pop * (srswor(n[1], N) > 0)
    j <- j + 1
  }
  for (i in 2:T) {
    unif <- runif(N, 0, 1)
    pop <- (unif <= lambda) * unif + pop
    if (t[j] == i) {
      data[, j] <- pop * (srswor(n[j], N) > 0)
      j <- j + 1
      if (j > lent) break
    }
  }
  return(data)
}

#' Simulate animal movement identification data with model C
#'
#' @param Z Total population of the whole area
#' @param N Initial population with the study area
#' @param T total time
#' @param n Number of identification in each observation
#' @param t Time of each observation if t is a list or array, the number
#'   of observation if t is an integer.
#' @param lambda move out rate
#' @param mu move in rate
#' @param seed random seed
#'
#' @return matrix of Z rows and length(t) columns
#' @export
#'
#' @examples
LIR.simulate.C <- function(Z, N, T, n, t, lambda, mu, seed = NULL) {
  set.seed(seed)
  lent <- length(t)
  if (lent == 1) warning("For model C, t had better be a list")
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
  t <- sort(t)
  if (1 == t[j]) {
    data[sample(pop_in, n[1]), 1] <- 1
    j <- j + 1
  }
  for (i in 2:T) {
    unif_in <- runif(length(pop_in))
    unif_out <- runif(length(pop_out))
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
    if (i == t[j]) {
      data[sample(pop_in, n[j]), j] <- 1
      j <- j + 1
      if (j > lent) break
    }
  }
  return(data)
}

#' Calculate lagged identification of each observation pair.
#' Number of observation is guessed by the column number of input matrix
#'
#' @param data Matrix of identification with row number equal to population
#'   and column number equal to observation.
#' @param t List like time of observation. Default NULL, must given if tau is true
#' @param tau whether return tau
#' @param n whether return n
#'
#' @return List of lagged identification number for each observation pair.
#' @export
#'
#' @examples
LIR.pairwise <- function(data, t = NULL, tau = TRUE, n = TRUE) {
  lent <- length(t)
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
  if (tau) {
    if (is.null(t))
      stop("No observation time provided to calculate tau")
    obs[["tau"]] <- unlist(
      sapply(1:(lent-1), FUN=function(i) {
        sapply((i+1):lent, FUN = function(j){return(t[j] - t[i])})
      })
    )
  }
  if (n) {
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

#' Perform a single bootstrap on the given observation matrix.
#' A SRSWR will be performed on individuals (row of the matrix)
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

