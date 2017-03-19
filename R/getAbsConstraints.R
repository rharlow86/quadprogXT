getAbsConstraints <- function(M, AmatL1){  
    k <- nrow(M)
    n <- ncol(M)
    nvar <- n  + k *2
    
    AEQ <- matrix(0, k, nvar)
    AEQ[ ,1:n] <- M
    AEQ[ ,-(1:n)] <- cbind(-diag(k), diag(k))

    AINEQ <- matrix(0, 2 * k, nvar)
    AINEQ[ ,-(1:n)] <- diag(2 * k)

    AINEQL1 <- matrix(0, ncol(AmatL1), nvar)
    AmatL1t <- t(AmatL1)
    AINEQL1[ ,-(1:n)] <- cbind(AmatL1t, AmatL1t)

    list(AEQ = AEQ, AINEQ = rbind(AINEQ, AINEQL1))
}

getAbsChangeConstraints <- function(M, AmatL1Change){
    stop("no function")
}
