nullConstraint <- function(){
    structure(list(Amat = NULL, bvec = NULL, meq = 0), class = "nullConstraint")
}

is.nullConstraint <- function(x){
    if(inherits(x, "nullConstraint")){
        return(TRUE)
    }else{
        return(FALSE)
    }
}

originalConstraints <- function(Amat, bvec, meq){
    if(is.null(Amat) || is.null(bvec) || is.null(meq)){
        nullConstraint()
    }else{
        structure(list(Amat = t(Amat), bvec = bvec, meq = meq), class = "originalConstraints")
    }
}

absConstraints <- function(M, AmatL1, bvecL1){

    if(!is.null(AmatL1) & !is.null(bvecL1)){

        k <- nrow(M)
        n <- ncol(M)
        nvar <- n + k * 2

        if(nrow(AmatL1) != k || length(bvecL1) != ncol(AmatL1)){
            stop("M, AmatL1 and bvecL1 are incompatible!")
        }

        AEQ <- matrix(0, k, nvar)
        AEQ[ ,1:n] <- M
        AEQ[ ,-(1:n)] <- cbind(-diag(k), diag(k))
        
        AINEQ <- matrix(0, 2 * k, nvar)
        AINEQ[ ,-(1:n)] <- diag(2 * k)
        
        AINEQL1 <- matrix(0, ncol(AmatL1), nvar)
        AmatL1t <- t(AmatL1)
        AINEQL1[ ,-(1:n)] <- cbind(AmatL1t, AmatL1t)
        
        Amat <- rbind(AEQ, AINEQ, AINEQL1)
        bvec <- c(rep(0,  k * 3), bvecL1)

        cons <- structure(
            list(Amat = Amat, bvec = bvec, meq = nrow(AEQ)),
            class = "absConstraints"
        )
        
        return(cons)
        
    }else{
        
        return(nullConstraint())
        
    }
           
}

absDeltaConstraints <- function(M, AmatL1Delta, bvecL1Delta, cvec, b0){

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
        
    }
    
    nullConstraint()
}
