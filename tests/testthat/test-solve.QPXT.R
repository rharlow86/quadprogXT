library(quadprogXT)

N <- 10
set.seed(2)
cr <- matrix(runif(N * N, 0, .05), N, N)
diag(cr) <- 1
cr <- (cr + t(cr)) / 2
set.seed(3)
sigs <- runif(N, min = .02, max = .25)
set.seed(5)

dvec <- runif(N, -.1, .1)
Dmat <- sigs %o% sigs * cr
Amat <- cbind(diag(N), diag(N) * -1)
bvec <- c(rep(-1, N), rep(-1, N))

#solve.QPXT(Dmat, dvec, Amat, bvec)

resBase <- solve.QPXT(Dmat, dvec, Amat, bvec)

test_that("call with factorized fails", {
    expect_error(solve.QPXT(Dmat, dvec, Amat, bvec, factorized = TRUE))
})

test_that("QPXT returns expected results for sum of absolute values <= 1 example", {
    res <- solve.QPXT(Dmat, dvec, Amat, bvec, AmatPosNeg = matrix(rep(-1, 2 * N)), bvecPosNeg = -1)
    expect_true(sum(abs(res$solution[1:N])) <= 1 + 1e-10)
})

test_that("QPXT still handles case where dvecPosNeg is not null (L1 norm penalty)", {
    resL1Penalty <- solve.QPXT(Dmat, dvec, Amat, bvec, dvecPosNeg = -.005 * rep(1, 2 * N))
    expect_true(sum(abs(resL1Penalty$solution[1:N]))  < sum(abs(resBase$solution)))
})

test_that("QPXT still handles case where dvecPosNeg is not null (L1 norm penalty)", {
    b0 <- rep(.15, N)
    thresh <- .25
    res <- solve.QPXT(Dmat, dvec, Amat, bvec, b0 = b0,
                      AmatPosNegDelta = matrix(rep(-1, 2 * N)), bvecPosNegDelta = -thresh)
    expect_true(sum(abs(res$solution[1:N] - b0)) <= thresh + 1e-10)
})
