expandConstraints <- function(cons, NABS, NABSCHANGE) UseMethod("expandConstraints")

#' @export
expandConstraints.originalConstraints <- function(cons, NABS, NABSCHANGE){
    cons$Amat <- cbind(cons$Amat, matrix(0, nrow(cons$Amat), NABS + NABSCHANGE))
    return(cons)
}

#' @export
expandConstraints.posNegConstraints <- function(cons, NABS, NABSCHANGE){
    cons$Amat <- cbind(cons$Amat, matrix(0, nrow(cons$Amat), NABSCHANGE))
    return(cons)
}

#' @export
expandConstraints.nullConstraint <- function(cons, NABS, NABSCHANGE){
    return(cons)
}
