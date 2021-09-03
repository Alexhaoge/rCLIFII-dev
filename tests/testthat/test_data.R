
test_that("Simulation A has correct shape and non-zero numbers", {
  n <- rep(c(3,4,5,1,2), 5)
  N <- 100
  tp <- c(1:5, 51:55, 101:105, 501:505, 601:605)
  data <- LIR.simulate.A(N, n, tp, seed = 0)
  expect_equal(dim(data), c(N, length(tp)))
  expect_equal(colSums(data != 0), n)
})

test_that("Simulation B has correct shape and non-zero numbers", {
  n <- rep(c(3,4,5,1,2), 5)
  N <- 100
  tp <- c(1:5, 51:55, 101:105, 501:505, 601:605)
  data <- LIR.simulate.B(N, n, tp, lambda=0.08)
  expect_equal(dim(data), c(N, length(tp)))
  expect_equal(colSums(data != 0), n)
  expect_error(LIR.simulate.B(N, c(n,n), tp, 0.01))
})

test_that("Simulation C", {
  # test length
  n <- rep(c(3,4,5,1,2), 5)
  N <- 100
  Z <- 300
  tp <- c(1:5, 51:55, 101:105, 501:505, 601:605)
  data <- LIR.simulate.C(Z, N, n, tp, lambda=0.08, mu=0.04, seed = 0)
  expect_equal(dim(data), c(Z, length(tp)))
  expect_equal(colSums(data != 0), n)
  expect_warning(LIR.simulate.C(Z, N, n, tp, lambda=0.1, mu=0.04))
  # test correctness by seed
  data <- matrix(0, 3, 3)
  data[1,1] = data[2,2] = data[1,3] = 1
  expect_equal(LIR.simulate.C(3, 1, 1, c(1,2,4), 0.08, 0.04, 1e9+7), data)
})

load('../data/ds.rda')

test_that("Pairwise and Non-parametric LIR", {
  # non.lir <- function(data, tp) {
  #   n <- colSums(data)
  #   len <- length(tp)
  #   tT <- max(tp)
  #   R.m <- rep(0, tT)
  #   R.n <- rep(0, tT)
  #   mij <- c()
  #   ni <- c()
  #   nj <- c()
  #   tauij <- c()
  #   k <- 1
  #   for (i in 1:(len-1)){
  #     for (j in (i+1):len){
  #       mij[k] <- sum((data[, i] == data[, j]) * (data[, i] != 0))
  #       ni[k] <- n[i]
  #       nj[k] <- n[j]
  #       tauij[k] <- tp[j] - tp[i]
  #       R.m[tauij[k]] <- R.m[tauij[k]] + mij[k]
  #       R.n[tauij[k]] <- n[i]*n[j] + R.n[tauij[k]]
  #       k <- k + 1
  #     }
  #   }
  #   R.tauij <- R.m[tauij]/R.n[tauij]
  #   tau <- sort(unique(tauij))
  #   R.tau <- R.m[tau]/R.n[tau]
  #   R.dat <- list(Rtau=Rtau, tau=tau, R.m=R.m, R.n=R.n,
  #                 m=mij, ni=ni, nj=nj, tauij=tauij)
  #   return(R.dat)
  # }
  tp <- c(1:3, 51:53)
  obs <- LIR.pairwise(ds, tp)
  load('../data/Rdat.rda')
  expect_equal(obs$m, Rdat$mij)
  expect_equal(obs$ni, Rdat$ni)
  expect_equal(obs$nj, Rdat$nj)
  expect_equal(obs$tauij, Rdat$tauij)
  expect_equal(obs$Rtau, Rdat$Rtau)
  expect_equal(obs$tau, Rdat$tau)
  rm(Rdat)
})

test_that("LIR bootstrap", {
  seed <- 1810081
  n <- nrow(ds)
  set.seed(seed)
  boot_row <- sample(1:n, n, replace = TRUE)
  boot <- LIR.bootstrap(ds, seed)
  expect_equal(boot, ds[boot_row, ])
})
