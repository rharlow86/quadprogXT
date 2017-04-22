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

res <- solve.QPXT(Dmat, dvec, Amat, bvec, AmatPosNeg = matrix(rep(-1, 2 * N)), bvecPosNeg = -1)


test_that("QPXT returns expected results for sum of absolute values <= 1 example", {
    expect_true(sum(abs(res$solution[1:N])) <= 1 + 1e-10)
})


