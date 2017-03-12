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
#' Amat * b >= bvec \cr
#' AmatL1 * |M * b| >= bvecL1 \cr
#' AmatL1Delta * |M * (b - b0)| >= bvecL1Delta \cr
#' 
#' 
#' @inheritParams quadprog::solve.QP
#'
#' @param M a k x n matrix that transforms the decision variable into a new space
#' to evaluate constraints. For most problems, this should just be the diagonal matrix.
#' 
#' @param AmatL1 m x k matrix of absolute value constraints on a linear mapping of b.
#' @param bvecL1 m x 1 vector of thresholds corresponding to AmatL1
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
                       b0 = NULL,
                       AmatL1Delta = NULL,
                       bvecL1Delta = NULL,
                       cvec = NULL
                       ){
    2
}
