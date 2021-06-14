#' Composite likelihood function for LIR(Lagged Identification Rate)
#'
#' @param theta Parameter to calculate R_tau
#' @param model Model to calculate R_tau
#' @param n number of individuals at each observation
#' @param N population
#'
#' @return Likelihood score(numeric). Inf if theta out of [0,1]
#' @export
#'
#' @examples
move.CL <- function(theta, model, n, N, m) {
  R_tau = model(theta)
  if (min(R_tau * n / N)<=0) return(Inf);
  if (max(R_tau * n / N)>=1) return(Inf);
  object <- sum(m * log(R_tau*n/N)+(n-m)*log(1-R_tau*n/N))
  return(-object)
}
