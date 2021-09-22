theta <- c(0.1, 0.1, 0.1)

test_that("test export model ABC", {
  expect_equal(LIR.model.A(theta[1], 0), theta[1])
  expect_equal(LIR.model.B(theta[1:2], 51), 6.0967e-4, tolerance=0.0001)
})

load('../data/ds.rda')
load('../data/Rdat.rda')
tp <- c(1:3, 51:53)

test_that("CL correctness", {
  cl_ans <- -27.8739108
  # expect_equal(LIR.CL(theta[1],'A',ds,tp,NULL), LIR.CL(theta[1],LIR.model.A,ds,tp,list(),-2.5))
  expect_equal(LIR.CL(theta[1],'A',ds,tp,NULL), cl_ans, tolerance=0.0001)
  expect_equal(LIR.CL(theta[1],LIR.model.A,ds,tp,list(),-2.5),
               cl_ans, tolerance=0.0001)
  expect_equal(LIR.CL(theta[1],'a',ds,tp,list(),Inf),
               cl_ans, tolerance=0.0001)
  expect_error(LIR.CL(theta[1],'xqwdsa', ds,tp, NULL))
  expect_equal(LIR.CL.pair(theta[1],LIR.model.A,Rdat$ni,
                           Rdat$nj,Rdat$mij,Rdat$tauij,NULL),
               cl_ans, tolerance=0.0001)

})

rm(ds)
rm(Rdat)
rm(tp)
