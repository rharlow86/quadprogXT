expandConstraints <- function(cons, NABS, NABSCHANGE) UseMethod("expandConstraints")

#' @export
expandConstraints.originalConstraints <- function(cons, NABS, NABSCHANGE){
    cons$Amat <- cbind(cons$Amat, matrix(0, nrow(cons$Amat), NABS + NABSCHANGE))
    return(cons)
}

#' @export
expandConstraints.L1 <- function(cons, NABS, NABSCHANGE){
    cons$Amat <- cbind(cons$Amat, matrix(0, nrow(cons$Amat), NABSCHANGE))
    return(cons)
}

#' @export
expandConstraints.L1Delta <- function(cons, NABS, NABSCHANGE){
    ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@5@"]]));##:ess-bp-end:##
    cons$Amat <- cbind(cons$Amat, matrix(0, nrow(cons$Amat), NABSCHANGE))
    return(cons)
}


#' @export
expandConstraints.nullConstraint <- function(cons, NABS, NABSCHANGE){
    return(cons)
}
