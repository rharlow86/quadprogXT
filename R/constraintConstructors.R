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

L1Constraints <- function(N, Amat, bvec, dvec, type, b0 = NULL){

    if(!is.null(Amat) || !is.null(dvec)){

        dummyCons <- dummyConstraints(N, b0)
        cons <- merge(dummyCons, N = N, Amat = Amat, bvec = bvec)
        class(cons) <- type
        return(cons)
        
    }else{
        
        return(nullConstraint())
        
    }

}

dummyConstraints <- function(N, b0 = NULL){
    nvar <- 3 * N
    diagN <- diag(N)
    
    AEQ <- matrix(0, N, nvar)
    AEQ[1:N, 1:N] <- diagN
    AEQ[ ,-(1:N)] <- cbind(-diagN, diagN)
    
    AINEQ <- matrix(0, 2 * N, nvar)
    AINEQ[ ,-(1:N)] <- diag(2 * N)

    if(is.null(b0)){
        b0 <- rep(0, N)
    }
    
    structure(
        list(Amat = rbind(AEQ, AINEQ), bvec = c(b0, rep(0, 2 * N)), meq = N),
        class = "dummyConstraints"
    )
}

#' @export
merge.dummyConstraints <- function(x, y = NULL, N, Amat = NULL, bvec = NULL, ...){
        
    if(!is.null(Amat)){
        
        if(ncol(Amat) != length(bvec)){
            stop("number of constraints not equal to number of thresholds")
        }
        
        Amat <- cbind(matrix(0, ncol(Amat), N), t(Amat))
        x$Amat <- rbind(x$Amat, Amat)
        x$bvec <- c(x$bvec, bvec)
    }
    
    return(x)
    
}
