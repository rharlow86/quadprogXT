findProblemType <- function(dvecPosNeg, AmatPosNeg, dvecPosNegDelta, AmatPosNegDelta){
    
    posNegCheck <- FALSE    
    if(!is.null(dvecPosNeg) || !is.null(AmatPosNeg)){
        posNegCheck <- TRUE
    }

    posNegDeltaCheck <- FALSE
    if(!is.null(dvecPosNegDelta) || !is.null(AmatPosNegDelta)){
        posNegDeltaCheck <- TRUE
    }

    
    if(posNegCheck && !posNegDeltaCheck){
        ans <- "posNegOnly"
    }else if(!posNegCheck && posNegDeltaCheck){
        ans <- "posNegDeltaOnly"
    }else if(posNegCheck && posNegDeltaCheck){
        ans <- "full"
    }else{
        ans <- "standardQP"
    }
    
    return(ans)
}
