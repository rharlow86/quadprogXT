singleDummyMat <- function(N, Amat, bvec, meq, absAmat, absBvec, b0 = NULL){
    nvar <- 3 * N
    diagN <- diag(N)

    ## dummy equalities/inequalities
    AEQ <- matrix(0, nvar, N)
    AEQ[1:N, 1:N] <- diagN
    AEQ[-(1:N),] <- rbind(-diagN, diagN)
    
    if(is.null(b0)){
        b0 <- rep(0, N)
    }
    beq <- b0

    AINEQ <- matrix(0, nvar, 2 * N)
    AINEQ[-(1:N), ] <- diag(2 * N)
    bineq <- rep(0, 2 * N)

    ## check/add existing constraints
    ncons <- length(bvec)
    if(ncons > 0){
        ineqIndex <- which(1:ncons > meq)
        eqIndex <- which(1:ncons <= meq)
        Amat <- rbind(
            Amat,
            matrix(0, 2 * N, ncol(Amat))
        )

        if(length(eqIndex)){        
            AEQ <- cbind(Amat[, eqIndex], AEQ)
            beq <- c(bvec[eqIndex], beq)
        }

        if(length(ineqIndex)){
            AINEQ <- cbind(Amat[ ,ineqIndex], AINEQ)
            bineq <- c(bvec[ineqIndex], bineq)
        }
        
    }

    ##check/add absolute value constraints
    if(!is.null(absAmat)){
        absAmat <- rbind(
            matrix(0, N, ncol(absAmat)),
            absAmat
        )

        AINEQ <- cbind(AINEQ, absAmat)
        bineq <- c(bineq, absBvec)
    }
    
    Amat <- cbind(AEQ, AINEQ)
    bvec <- c(beq, bineq)
    meq <- ncol(AEQ)
    
    list(Amat = Amat, bvec = bvec, meq = meq)
    
}
