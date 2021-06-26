
test_that("Simulation A has correct shape and non-zero numbers", {
  n <- rep(c(3,4,5,1,2), 5)
  N <- 100
  T <- 605
  t <- c(1:5, 51:55, 101:105, 501:505, 601:605)
  data <- LIR.simulate.A(N, n, t, seed = 0)
  expect_equal(dim(data), c(N, length(t)))
  expect_equal(colSums(data != 0), n)
})

test_that("Simulation B has correct shape and non-zero numbers", {
  n <- rep(c(3,4,5,1,2), 5)
  N <- 100
  T <- 605
  t <- c(1:5, 51:55, 101:105, 501:505, 601:605)
  data <- LIR.simulate.B(N, T, n, t, lambda=0.08)
  expect_equal(dim(data), c(N, length(t)))
  expect_equal(colSums(data != 0), n)
  expect_error(LIR.simulate.B(N, T, c(n,n), t, 0.01))
})

test_that("Simulation C", {
  n <- rep(c(3,4,5,1,2), 5)
  N <- 100
  Z <- 300
  T <- 605
  t <- c(1:5, 51:55, 101:105, 501:505, 601:605)
  data <- LIR.simulate.C(Z, N, T, n, t, lambda=0.08, mu=0.04, seed = 0)
  expect_equal(dim(data), c(Z, length(t)))
  expect_equal(colSums(data != 0), n)
  expect_warning(LIR.simulate.C(Z, N, T, n, t, lambda=0.1, mu=0.04, seed = 0))
})

test_that("Pairwise", {
  t <- c(1:5, 51:55, 101:105, 501:505, 601:605)
  tau <- c()
  k <- 1
  for (i in 1:(length(t)-1)){
    for (j in (i+1):length(t)){
      tau[k] <- t[j] - t[i]
      k <- k + 1
    }
  }
  data <- replicate(25, c(0,1), simplify = 'matrix')
  m <- rep(1, 25 * 24 /2)
  obs <- LIR.pairwise(data, t)
  expect_equal(obs$m, m)
  expect_equal(obs$ni, m)
  expect_equal(obs$nj, m)
  expect_equal(obs$tau, tau)
  expect_error(LIR.pairwise(data, tau=T))
  obs <- LIR.pairwise(data, t, tau = F, n = F)
  expect_null(obs$tau)
  expect_null(obs$ni)
})
