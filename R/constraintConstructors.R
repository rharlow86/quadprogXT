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

posNegConstraints <- function(AmatPosNeg, bvecPosNeg){

    if(!is.null(AmatPosNeg) & !is.null(bvecPosNeg)){

        N <- nrow(AmatPosNeg) / 2
        nvar <- 3 * N

        if(ncol(AmatPosNeg) != length(bvecPosNeg)){
            stop("AmatPosNeg and bvecPosNeg are incompatible!")
        }

        AEQ <- matrix(0, N, nvar)
        AEQ[1:N, 1:N] <- diag(N)
        AEQ[ ,-(1:N)] <- cbind(-diag(N), diag(N))
        
        AINEQ <- matrix(0, 2 * N, nvar)
        AINEQ[ ,-(1:N)] <- diag(2 * N)
        
        AINEQPosNeg <- t(AmatPosNeg)
        AINEQPosNeg <- cbind(matrix(0, ncol(AmatPosNeg), N), AINEQPosNeg)
        
        
        Amat <- rbind(AEQ, AINEQ, AINEQPosNeg)
        bvec <- c(rep(0,  nvar), bvecPosNeg)

        cons <- structure(
            list(Amat = Amat, bvec = bvec, meq = nrow(AEQ)),
            class = "posNegConstraints"
        )
        
        return(cons)
        
    }else{
        
        return(nullConstraint())
        
    }
           
}
