merge.quadprogXTConstraintList <- function(x, y = NULL, N, ...){
    
    consClasses <- sapply(x, class)
    NABS <- any("L1" %in% consClasses) * 2 * N
    NABSCHANGE <- any("L1Delta" %in% consClasses) * 2 * N

    constraintListExp <- lapply(x, expandConstraints, N, NABS, NABSCHANGE)
    equalityConstraints <- lapply(constraintListExp, equalityConstraints)
    inequalityConstraints <- lapply(constraintListExp, inequalityConstraints)
    
    AEQ <- do.call("rbind", lapply(equalityConstraints, getElement, "Amat"))
    BEQ <- do.call("c", lapply(equalityConstraints, getElement, "bvec"))

    AINEQ <- do.call("rbind", lapply(inequalityConstraints, getElement, "Amat"))
    BINEQ <- do.call("c", lapply(inequalityConstraints, getElement, "bvec"))

    list(Amat = t(rbind(AEQ, AINEQ)), bvec = c(BEQ, BINEQ), meq = length(BEQ))
    
}
