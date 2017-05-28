#' Solve a quadratic program with absolute values in constraints & objective
#'
#' @description This function allows for absolute value constraints and absolute values in the
#' objective.  Additionally, this function implicitly takes advantage of sparsity in the constraint
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

#' @details In order to handle constraints on b_positive and b_negative, slack variables are introduced.  The total number of parameters in the problem increases by the following amounts: \cr
#' If all the new parameters (those not already used by quadprog) remain NULL, the problem size does not increase and quadprog::solve.QP is called after normalizing the constraint matrix and converting to a sparse matrix representation.\cr
#' If AmatPosNeg, bvecPosNeg or dvecPosNeg are not null, the problem size increases by 2 * n
#'
#' Despite the potential large increases in problem size, the underlying solver is written in
#' Fortran and converges quickly for problems involving even hundreds of parameters.  Additionally,
#' it has been the author's experience that solutions solved via the convex quadprog are much more
#' stable than those solved by other methods (e.g. a non-linear solver).
#'
#' @usage solve.QPXT(Dmat, dvec, Amat, bvec, meq = 0, factorized = FALSE)
#' @export solve.QPXT
solve.QPXT <- function(Dmat, dvec, Amat, bvec, meq = 0, factorized = FALSE,
                       AmatPosNeg = NULL,
                       bvecPosNeg = NULL,
                       dvecPosNeg = NULL,
                       b0 = NULL,
                       AmatPosNegDelta = NULL,
                       bvecPosNegDelta = NULL,
                       dvecPosNegDelta = NULL
                       ){

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
    diag(DMAT) <- 1e-8

    if(NVAR > N){
        if(is.null(dvecPosNeg)){
            dvecPosNeg <- rep(0, NVAR - N)
        }        
    }

    DMAT[1:N, 1:N] <- Dmat
    DVEC <- c(dvec, dvecPosNeg)

    comp <- convertToCompact(Amat)    
    res <- quadprog::solve.QP.compact(DMAT, DVEC, comp$Amat, comp$Aind, bvec, meq)

    return(res)
}
