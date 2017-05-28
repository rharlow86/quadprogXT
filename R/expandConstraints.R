expandConstraints <- function(cons, N, NABS, NABSCHANGE) UseMethod("expandConstraints")

#' @export
expandConstraints.originalConstraints <- function(cons, N, NABS, NABSCHANGE){
    cons$Amat <- cbind(cons$Amat, matrix(0, nrow(cons$Amat), NABS + NABSCHANGE))
    return(cons)
}

#' @export
expandConstraints.L1 <- function(cons, N, NABS, NABSCHANGE){
    cons$Amat <- cbind(cons$Amat, matrix(0, nrow(cons$Amat), NABSCHANGE))
    return(cons)
}

#' @export
expandConstraints.L1Delta <- function(cons, N, NABS, NABSCHANGE){    
    cons$Amat <- with(cons, cbind(Amat[ , 1:N], matrix(0, nrow(Amat), NABS), Amat[ ,-(1:N)]))
    return(cons)
}


#' @export
expandConstraints.nullConstraint <- function(cons, N, NABS, NABSCHANGE){
    return(cons)
}
