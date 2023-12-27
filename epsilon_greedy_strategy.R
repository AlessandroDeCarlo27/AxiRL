epsilon_greedy_strategy <- function(df_i,states_action_i,q_i,eps){
    
    actions <- states_action_i[states_action_i$ID==df_i$state_id,-(1:5)]
    r=runif(1)
    if (r<=eps) {
        #### random choice
        actions <- as.matrix(actions[,as.numeric(actions[1,])>0])
        idx_actions <- sample(1:ncol(actions),size = 1,prob = 1/actions) 
    }else{
        vals <- q_i[df_i$state_id,]
        names(vals) <- colnames(actions)
        vals<- vals[as.numeric(actions[1,])>0]
        actions <- as.matrix(actions[,as.numeric(actions[1,])>0])
        maxV <- max(vals)
        idx <- which(vals==maxV)
        if (length(idx)==1) {
            idx_actions <- idx[1]
        }else{
            idx_actions <- sample(1:ncol(actions),size = 1,prob = 1/actions)   
        }
    }
    
    name_action <- colnames(actions)[idx_actions]
    return(name_action)
}