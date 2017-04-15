merge.quadprogXTConstraintList <- function(x, y = NULL, N, K, ...){
    
    consClasses <- sapply(x, class)
    NABS <- any("absConstraints" %in% consClasses) * 2 * K
    NABSCHANGE <- any("absDeltaConstraints" %in% consClasses) * 2 * K

    constraintListExp <- lapply(x, expandConstraints, NABS, NABSCHANGE)
    equalityConstraints <- lapply(constraintListExp, equalityConstraints)
    inequalityConstraints <- lapply(constraintListExp, inequalityConstraints)
    
    AEQ <- do.call("rbind", lapply(equalityConstraints, getElement, "Amat"))
    BEQ <- do.call("c", lapply(equalityConstraints, getElement, "bvec"))

    AINEQ <- do.call("rbind", lapply(inequalityConstraints, getElement, "Amat"))
    BINEQ <- do.call("c", lapply(inequalityConstraints, getElement, "bvec"))

    list(Amat = t(rbind(AEQ, AINEQ)), bvec = c(BEQ, BINEQ), meq = length(BEQ))
    
}
