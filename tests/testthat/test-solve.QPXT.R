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

res <- solve.QPXT(Dmat, dvec, Amat, bvec, M = M, AmatL1 = AmatL1, bvecL1 = -1)

eres <- c(1.05390696170105e-12, -5.17394481433108e-13, 0.130657802509203, 
-2.10252600911269e-14, -0.439638172114367, -5.49495825072432e-13, 
-9.13653837593511e-14, -1.41453494329112e-12, 0.0146144040327969, 
-0.415089621337841, -2.35334549462209e-13, 1.66761364963982e-14, 
0.13065780250948, 3.65655179749206e-13, -3.70727291536499e-13, 
-9.52001931637194e-17, 2.93461327956727e-14, -6.43354410819208e-16, 
0.0146144040321515, -5.69286923172036e-13, -2.29522675021771e-16, 
0, 5.34048920231976e-13, 6.31419536253527e-13, 1.69629210507801e-16, 
0.439638172114357, -3.73376449323299e-13, 6.96666898844059e-14, 
1.01005678432769e-12, -1.31295827364576e-13, 0.415089621337312, 
4.51461724332059e-13)

test_that("QPXT returns expected results for AmatL1 example", {
    expect_equal(res$solution, eres)
})


