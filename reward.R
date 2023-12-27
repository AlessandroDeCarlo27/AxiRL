get_reward_dbp<- function(x){
    
    if (x<60) {
        return(0.1)
    }else if(x>=60 & x<90){
        k21 <- log(0.1/5)*(-1/30)
        y21 <- 5*exp(-k21*abs(x-90))
        return(y21)
    }else if(x>=90 & x<92.5){
        m2 <- 3/2.5
        y2 <- m2*x+(7-(m2*90))
        return(y2)
    }else if(x>=92.5 & x<97.5){
        return(10)
    }else if(x>=97.5 & x<100){
        m2 <- 3/2.5
        y4 <- -m2*x+(7+(m2*100))
        return(y4)
    }else if(x>=100 & x<=105){
        k5 <- -0.2*log(0.1/5)
        y5 = 5*exp(-k5*(x-100))
        return(y5)
    }else if(x>105){
        return(0.1)
    }
}

get_reward_AX <- function(bp,ax){
    
    if (bp<100 & ax>0) {
            doses<-c(2,3,5,7,10)
            x <- which(doses==ax) 
            return(x*2)
    }else{
        return(0)
    }
}

get_reward_AH <- function(bp,ah){
    
    if (bp<100 & bp>=90) {
        doses <- seq(4,0,-1)
        return(which(doses==ah))
    }else{
        return(0)
    }
}

get_reward <- function(df_i){
    
    rew <-get_reward_dbp(as.numeric(df_i["bp"])) + 
        get_reward_AX(as.numeric(df_i["bp"]),as.numeric(df_i["prevAX"])) +
        get_reward_AH(as.numeric(df_i["bp"]),as.numeric(df_i["prevAH"]))
    return(rew)
}