equalityConstraints <- function(cons){
    if(cons$meq > 0){
        idx <- 1:cons$meq
        return(list(Amat = cons$Amat[idx, ], bvec = cons$bvec[idx]))
    }else{
        return(list(Amat = NULL, bvec = NULL))
    }
}

inequalityConstraints <- function(cons){
    if(!is.nullConstraint(cons)){
        if(cons$meq > 0){
            idx <- 1:cons$meq
            return(list(Amat = cons$Amat[-idx, ], bvec = cons$bvec[-idx]))
        }else{
            return(cons)
        }
    }else{
        return(list(Amat = NULL, bvec = NULL))
    }
}

