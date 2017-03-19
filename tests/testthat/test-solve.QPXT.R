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

M <- diag(N)
M <- rbind(M, 0)
set.seed(11)
M[nrow(M), ] <- runif(N)
AmatL1 <- matrix(-1, nrow(M), 1)

dvecl1 <- rep(0, 2 * nrow(M))

res <- solve.QPXT(Dmat, dvec, Amat, bvec, M = M, AmatL1 = AmatL1, bvecL1 = -1, dvecL1 = dvecl1)

eres <- c(-2.57849297469193e-13, 0.817552559965783, 1, -3.190404458172e-13, 
-0.554422073014709, 1, 1, 1, 1, -0.445577926987129, -1.06438693014905e-13, 
1.05332496803522e-12, 2.7649140658094e-13, -8.95163249580714e-14, 
2.68581731555996e-14, -2.12877412522672e-13, -1.32717689010987e-13, 
2.59192423708128e-13, -7.25384867547912e-14, -1.27061561261557e-13, 
-1.52174624590743e-13, 0, 1.09349896511816e-13, 5.97856496093912e-13, 
0, 0.554422073014355, -6.96828998504559e-14, -2.73787766597224e-13, 
-6.79929736839417e-14, 1.2446894793319e-13, 0.445577926986847, 
2.13737765001684e-16)

test_that("QPXT returns expected results for AmatL1 example", {
    expect_equal(res$solution, eres)
})
