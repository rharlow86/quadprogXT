#' Solve a quadratic program with absolute values in constraints & objective
#'
#' @description This function allows for absolute values in the constraint matrix and in the
#' objective.  Additionally, this function implicitly takes advantage of sparsity in the constraint
#' matrix and can improve numerical stability by normalizing the constraint matrix. For the
#' rest of the documentation, assume that Dmat is n x n.\cr
#' \preformatted{
#' 
#' The solver solves the following problem (each * corresponds to matrix multiplication):
#'
#' min( -t(dvec) * b + 1/2 t(b) * Dmat * b + t(cvec) * |M * (b - b0)|)
#' }
#' \cr
#' s.t. \cr
#' t(Amat) * b >= bvec \cr
#' t(AmatL1) * |M * b| >= bvecL1 \cr
#' t(AmatL1Delta) * |M * (b - b0)| >= bvecL1Delta \cr
#' 
#' 
#' @inheritParams quadprog::solve.QP
#'
#' @param M a k x n matrix that transforms the decision variable into a new space
#' to evaluate constraints. For most problems, this should just be the diagonal matrix.
#' 
#' @param AmatL1 k x m matrix of absolute value constraints on a linear mapping of b.
#' @param bvecL1 m length vector of thresholds corresponding to AmatL1
#' 
#' @param b0 a starting point that describes the 'current' state of the problem such that
#' constraints and penalty on absolute changes in the decision variable from a starting point can
#' be incorporated.  b0 is an n x 1 vector. Note that b0 is NOT a starting point for the
#' optimization - that is handled implicitly by quadprog.
#' 
#' @param AmatL1Delta m2 x k matrix of absolute value constraints on a linear mapping of changes
#' in b.
#' @param bvecL1Delta m2 x 1 vector of thresholds corresponding to AmatL1Delta
#'
#' @param cvec a k x 1 vector of 'costs' associated with absolute changes in a mapping of changes
#' in b.
#' 
#' @details In order to handle absolute value constraints, slack variables are introduced.  The
#' total number of parameters in the problem increases by the following amounts: \cr
#' If all the new parameters remain NULL, the problem size does not increase and quadprog is called
#' after normalizing the constraint matrix and converting to a sparse matrix representation.\cr
#' If b0, AmatL1Delta, bvecL1Delta and cvec are all null, the problem increases in size by 2 * k.\cr
#' If all new parameters are not null, the problem size increases by 4 * k. \cr
#'
#' Despite the potential large increases in problem size, the underlying solver is written in
#' Fortran and converges quickly for problems involving even hundreds of parameters.  Additionally,
#' it has been the author's experience that solutions solved via the convex quadprog are much more
#' stable than those solved by other methods (e.g. a non-linear solver).
#'
#' @usage solve.QPXT(Dmat, dvec, Amat, bvec, meq = 0, factorized = FALSE)
#' @export solve.QPXT
solve.QPXT <- function(Dmat, dvec, Amat, bvec, meq = 0, factorized = FALSE,
                       M = NULL,
                       AmatL1 = NULL,
                       bvecL1 = NULL,
                       dvecL1 = NULL,
                       b0 = NULL,
                       AmatL1Delta = NULL,
                       bvecL1Delta = NULL,
                       cvec = NULL
                       ){

    N <- length(dvec)
    NDUMMYABS <- NDUMMYABSCHANGE <- NDUMMY <- K <- 0    

    if(!is.null(M)) {
        K <- nrow(M)
    }

    if(!is.null(AmatL1) && !is.null(bvecL1) && !is.null(dvecL1)){        

        if(nrow(AmatL1) != K || length(bvecL1) != ncol(AmatL1) || length(dvecL1) != (2 * K)){
            stop("M, AmatL1 and bvecL1 are incompatible!")
        }
        
        absConstraints <- getAbsConstraints(M, AmatL1)
        NDUMMYABS <- 2 * K
    }

    if(!is.null(cvec) || !is.null(AmatL1Delta)){

        if(!is.null(cvec)){
            
            if(length(cvec) != K || length(b0) != N){
                stop("cvec, M and b0 are incompatible")
            }
        }

        if(!is.null(AmatL1Delta)){

            if((ncol(AmatL1Delta) != length(bvecL1Delta)) || length(b0) != N){
                stop("AmatL1Delta, bvecL1Delta and b0 are incompatible")
            }
            
        }

        absChangeConstraints <- getAbsChangeConstraints(M, AmatL1Delta)
        NDUMMYABSCHANGE <- 2 * K        
    }

    NDUMMY <- sum(NDUMMYABS, NDUMMYABSCHANGE)
    NVAR <- N + NDUMMY
    
    if(NDUMMY){
        ##more used to rows to think of constraints as rows
        Amat <- t(Amat)
        NINEQ <- nrow(Amat) - meq
        BINEQ <- BEQ <- NULL
        
        AEQ <- matrix(0, meq, NVAR) 
        AINEQ <- matrix(0, NINEQ, NVAR)
        
        if(meq){
            AEQ[1:meq, 1:N] <- Amat[1:meq, 1:N]
            BEQ <- bvec[1:meq]
        }
        
        if(NINEQ){
            BINEQ <- bvec[(meq+1) : length(bvec)]
            AINEQ[1:NINEQ, 1:N] <- Amat[(meq + 1):nrow(Amat), ]
        }
                
        if(length(absConstraints)){
            AEQ <- rbind(AEQ, absConstraints$AEQ)
            AINEQ <- rbind(AINEQ, absConstraints$AINEQ)
            BEQ <- c(BEQ, rep(0, nrow(absConstraints$AEQ)))
            BINEQ <- c(BINEQ, rep(0, NDUMMYABS), bvecL1)            
        }
        
        Amat <- rbind(AEQ, AINEQ)
        bvec <- c(BEQ, BINEQ)
        Amat <- t(Amat)
        
    }

    
    norms <- normalizeConstraints(Amat, bvec)
    Amat <- norms$Amat
    bvec <- norms$bvec

    DMAT <- matrix(0, NVAR, NVAR)
    diag(DMAT) <- 1e-8

    DMAT[1:N, 1:N] <- Dmat
    DVEC <- c(dvec, dvecL1)

    comp <- convertToCompact(Amat)    
    res <- quadprog::solve.QP.compact(DMAT, DVEC, comp$Amat, comp$Aind, bvec, meq)

    return(res)
}
