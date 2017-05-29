#' Solve a quadratic program with absolute values in constraints & objective
#'
#' @description solveQPXT  allows for absolute value constraints and absolute values in the
#' objective.  buildQP builds a parameter list that can then be passed to
#' quadprog::solve.QP.compact directly if desired by the user.
#' solveQPXT implicitly takes advantage of sparsity in the constraint
#' matrix and can improve numerical stability by normalizing the constraint matrix. For the
#' rest of the documentation, assume that Dmat is n x n.\cr
#' 
#' The solver solves the following problem (each * corresponds to matrix multiplication):
#' \preformatted{
#' min(
#' -t(dvec) * b + 1/2 t(b) * Dmat * b + -t(dvecPosNeg) * c(b_positive, b_negative)
#'  + -t(dvecPosNegDelta) * c(deltab_positive, deltab_negative)
#' )
#' s.t.
#' t(Amat) * b >= bvec 
#' t(AmatPosNeg) * c(b_positive, b_negative) >= bvecPosNeg
#' t(AmatPosNegChange) * c(deltab_positive, deltab_negative) >= bvecPosNegChange
#' b_positive, b_negative >= 0, b = b_positive - b_negative
#' deltab_positive >= 0, deltab_negative >= 0, b - b0 = deltab_positive - deltab_negative
#' }
#' 
#' @inheritParams quadprog::solve.QP
#'
#' 
#' @param AmatPosNeg k x 2n matrix of constraints on the positive and negative part of b
#' @param bvecPosNeg k length vector of thresholds to the constraints in AmatPosNeg
#' @param dvecPosNeg k x 2n vector of loadings on the positive and negative part of b, respectively
#' 
#' @param b0 a starting point that describes the 'current' state of the problem such that
#' constraints and penalty on absolute changes in the decision variable from a starting point can
#' be incorporated.  b0 is an n x 1 vector. Note that b0 is NOT a starting point for the
#' optimization - that is handled implicitly by quadprog.
#' @param AmatPosNegDelta l x 2n matrix of constraints on the positive and negative part of a change in b from a starting point, b0.
#' @param bvecPosNegDelta l length vector of thresholds to the constraints in AmatPosNegDelta
#' @param dvecPosNegDelta l x 2n vector of loadings in the objective function on the positive and negative part of changes in b from a starting point of b0.
#' @param tol tolerance along the diagonal of the expanded Dmat for slack variables
#' @details In order to handle constraints on b_positive and b_negative, slack variables are introduced.  The total number of parameters in the problem increases by the following amounts: \cr
#' If all the new parameters (those not already used by quadprog) remain NULL, the problem size does not increase and quadprog::solve.QP is called after normalizing the constraint matrix and converting to a sparse matrix representation.\cr
#' If AmatPosNeg, bvecPosNeg or dvecPosNeg are not null, the problem size increases by 2 * n
#' If AmatPosNegDelta or devecPosNegDelta are not null, the problem size increases by 2 * n.
#' This results in a potential problem size of up to 5 * n.
#' Despite the potential large increases in problem size, the underlying solver is written in
#' Fortran and converges quickly for problems involving even hundreds of parameters.  Additionally,
#' it has been the author's experience that solutions solved via the convex quadprog are much more
#' stable than those solved by other methods (e.g. a non-linear solver).
#'
#' Note that due to the fact that the constraints are normalized, the original constraint values the user passed will not be returned by buildQP. 
#' @export solveQPXT
#'
#' @examples
#' ##quadprog example"
#' Dmat       <- matrix(0,3,3)
#' diag(Dmat) <- 1
#' dvec       <- c(0,5,0)
#' Amat       <- matrix(c(-4,-3,0,2,1,0,0,-2,1),3,3)
#' bvec       <- c(-8,2,0)
#' qp <- quadprog::solve.QP(Dmat,dvec,Amat,bvec=bvec)
#' qpXT <- solveQPXT(Dmat,dvec,Amat,bvec=bvec)
#' range(qp$solution - qpXT$solution)
#' 
#' N <- 10
#' set.seed(2)
#' cr <- matrix(runif(N * N, 0, .05), N, N)
#' diag(cr) <- 1
#' cr <- (cr + t(cr)) / 2
#' set.seed(3)
#' sigs <- runif(N, min = .02, max = .25)
#' set.seed(5)

#' dvec <- runif(N, -.1, .1)
#' Dmat <- sigs %o% sigs * cr
#' Amat <- cbind(diag(N), diag(N) * -1)
#' bvec <- c(rep(-1, N), rep(-1, N))

#' resBase <- solveQPXT(Dmat, dvec, Amat, bvec)
#' ##absolute value constraint on decision variable:
#' res <- solveQPXT(Dmat, dvec, Amat, bvec, AmatPosNeg = matrix(rep(-1, 2 * N)), bvecPosNeg = -1)
#' sum(abs(res$solution[1:N]))
#'
#' ## penalty of L1 norm
#' resL1Penalty <- solveQPXT(Dmat, dvec, Amat, bvec, dvecPosNeg = -.005 * rep(1, 2 * N))
#' sum(abs(resL1Penalty$solution[1:N]))
#'
#' ## constraint on amount decision variable can vary from a starting point
#' b0 <- rep(.15, N)
#' thresh <- .25
#' res <- solveQPXT(Dmat, dvec, Amat, bvec, b0 = b0,
#' AmatPosNegDelta = matrix(rep(-1, 2 * N)), bvecPosNegDelta = -thresh)
#' sum(abs(res$solution[1:N] - b0))
#'
#' ##use buildQP, then call solve.QP.compact directly
#' qp <- buildQP(Dmat, dvec, Amat, bvec, b0 = b0,
#' AmatPosNegDelta = matrix(rep(-1, 2 * N)), bvecPosNegDelta = -thresh)
#' res2 <- do.call(quadprog::solve.QP.compact, qp)
#' range(res$solution - res2$solution)

solveQPXT <- function(Dmat, dvec, Amat, bvec, meq = 0, factorized = FALSE,
                       AmatPosNeg = NULL,
                       bvecPosNeg = NULL,
                       dvecPosNeg = NULL,
                       b0 = NULL,
                       AmatPosNegDelta = NULL,
                       bvecPosNegDelta = NULL,
                       dvecPosNegDelta = NULL,
                       tol = 1e-8
                       ){
    
    args <- as.list(environment())
    qpArgs <- do.call(buildQP, args)
    res <- do.call(quadprog::solve.QP.compact, qpArgs)
    return(res)
}

#' @rdname solveQPXT
#' @export
buildQP <- function(Dmat, dvec, Amat, bvec, meq = 0, factorized = FALSE,
                    AmatPosNeg = NULL,
                    bvecPosNeg = NULL,
                    dvecPosNeg = NULL,
                    b0 = NULL,
                    AmatPosNegDelta = NULL,
                    bvecPosNegDelta = NULL,
                    dvecPosNegDelta = NULL,
                    tol = 1e-8
                    ){
    
    if(factorized){
        stop("solveQPXT does not handle a factorized Dmat (yet)")
    }
    
    N <- length(dvec)
    
    constraintList <- structure(
        list(
            originalConstraints(Amat, bvec, meq),
            L1Constraints(N, AmatPosNeg, bvecPosNeg, dvecPosNeg, "L1"),
            L1Constraints(N, AmatPosNegDelta, bvecPosNegDelta, dvecPosNegDelta, "L1Delta",b0)
        ),
        class = "quadprogXTConstraintList"
    )
    
    constraints <- merge(constraintList, N = N)    
    norms <- normalizeConstraints(constraints$Amat, constraints$bvec)

    Amat <- constraints$Amat
    bvec <- constraints$bvec
    meq <- constraints$meq
    
    NVAR <- nrow(Amat)
    DMAT <- matrix(0, NVAR, NVAR)
    diag(DMAT) <- tol

    if(NVAR > N){
        if(is.null(dvecPosNeg)){
            dvecPosNeg <- rep(0, 2 * N * !is.null(constraintList[[2]]$Amat))
        }

        if(is.null(dvecPosNegDelta)){
            dvecPosNegDelta <- rep(0, 2 * N * !is.null(constraintList[[3]]$Amat))
        }
    }

    DMAT[1:N, 1:N] <- Dmat
    DVEC <- c(dvec, dvecPosNeg, dvecPosNegDelta)

    comp <- convertToCompact(Amat)    

    list(Dmat = DMAT, dvec = DVEC, Amat = comp$Amat, Aind = comp$Aind, bvec = bvec, meq = meq,
         factorized = factorized)

}
