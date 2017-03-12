#' 'Sparsify' constraint matrix
#'
#' @param Amat a constraint matrix as defined in solve.QP
#' @return a list with two elements: Amat and Aind as necessary to be passed to solve.QP.compact
#' @seealso quadprog::solve.QP
#' @seealso quadprog::solve.QP.compact
#' @export
convertToCompact <- function(Amat){
    m_i <- colSums(abs(Amat) != 0)
    maxmi <- max(m_i)
    NC <- ncol(Amat)
    Aind <- matrix(0, maxmi + 1, NC)

    Aind[1, ] <- maxmi
    nonZero <- which(Amat != 0, arr.ind = TRUE)
    Aind[-1, ] <- nonZero[ , "row"]

    AmatSparse <- matrix(0, maxmi, NC)    
    for(i in 1:NC){
        iNonZero <- which(Amat[,i] != 0)
        if(!length(iNonZero)){
            stop(sprintf("There are not any non zero elements in column %i", i))
        }
        AmatSparse[ ,i] <- Amat[iNonZero ,i]
    }
    list(Amat = AmatSparse, Aind = Aind)
}

