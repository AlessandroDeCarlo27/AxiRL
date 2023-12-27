get_discrDBP <- function(x,isInitial){
    
    if (isInitial) {
        return(initDiscr(x))
    }else{
        return(treatDiscr(x))
    }
}

initDiscr <- function(x){
    if(x<70){
        return(-3)
    }else if(x>=70 & x<80){
        return(-2)
    }else{
        return(-1)
    }
}


treatDiscr <- function(x){
    
    if (x<90) {
        return(1)
    }else if (x>=90 & x<100) {
        return(2)
    }else if (x>=100 & x<=105){
        return(3)
    }else if (x>105 & x<=120){
        return(4)
    }else if(x>120){
        return(5)
    }
}

get_state_id<- function(x,states_map){
    
    id<- states_map$ID[(states_map$BPLevel==as.numeric(x["discr_bp"]) &
                   states_map$AX == as.numeric(x["prevAX"]) &
                   states_map$AH == as.numeric(x["prevAH"]) &
                   states_map$FlagInt == as.numeric(x["FlagInt"]))]
    return(id)
}

parse_dosing <- function(x){
    info <- unlist(strsplit(x,"_"))
    ax <- substr(info[1],3,nchar(info[1]))
    ah <- substr(info[2],3,nchar(info[2]))
    o <- list()
    o$AX <- ax
    o$AH <- ah
    return(o)
    
}


dosing_2_id <- function(x,states_map){
    
    cn <- colnames(states_action)[-(1:5)]
    return(which(cn==x))
}