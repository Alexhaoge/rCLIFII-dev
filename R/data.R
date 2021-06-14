
#' Simulate animal movement identification data with model A
#' Population is constant and no migration.
#'
#' @param N population
#' @param n number of observed individual at each moment
#' @param t number of observation
#' @param seed random seed
#'
#' @return 0-1 identification matrix
#' @export
#'
#' @examples
move.simulate.A <- function(N, n, t, seed = NULL) {
  set.seed(seed)
  data <- replicate(t, srswor(n, N))
  return(data)
}


move.pairwise.A <- function(data) {
  t = ncol(data)
  obs <- data.frame(
    m = unlist(
      sapply(1:(t-1), FUN=function(i) {
        sapply((i+1):t, FUN = function(j){sum(data[,i]*data[,j])})
      })
    )
    tau = unlist(
      sapply(1:(t-1), FUN=function(i) {
        sapply((i+1):t, FUN = function(j){j-i})
      })
    )
  )
  return(obs)
}

move.simulate.B <- function(N, n, t, seed = NULL) {
  set.seed(seed)

}
