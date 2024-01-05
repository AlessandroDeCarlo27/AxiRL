source("utils.R")
source("epsilon_greedy_strategy.R")
source("reward.R")

test_QL_agents <- function(patient_state_initial,patients,eps,states_action_list,
                            Q_table_ind.list,n_cyc){
    
    
    
    bp_list <- list()
    dose_list <- list()
    
    Npat <- nrow(patient_state_initial)
    discRew <- rep(0,Npat)
    patient_state <- patient_state_initial #set initial state
    patient_df <- patients
    
    for (i in 1:Npat) {
        df_i<-patient_state[i,]
        name_action <- epsilon_greedy_strategy(df_i=df_i,
                                               states_action_i = states_action_list[[i]],
                                               q_i = Q_table_ind.list[[i]],
                                               eps=eps)
        
        dose_info <- parse_dosing(name_action)
        patient_df$AX[i] <- dose_info$AX
        patient_df$AH[i] <- dose_info$AH
        if (eps>0) {
            states_action_list[[i]][states_action_list[[i]]$ID==df_i$state_id,name_action] <-
                states_action_list[[i]][states_action_list[[i]]$ID==df_i$state_id,name_action]+1    
        }
        patient_state$CurrentAct[i] <- name_action
    }
    
    
    for (j in 1:n_cyc) {
        
        o <- list(name=c('bp'),time=seq(0,28,by=1))
        res<- simulx(model = "dbp_model.txt",parameter = patient_df, output = o)
        # store output
        bp_tmp <- res$bp
        
        if (j>1) {
            bp_st <- res$bp %>% filter(time>0) %>% mutate(time=time+((j-1)*28))
        }else{
            bp_st <- bp_tmp
        }
        
        bp_list[[length(bp_list)+1]] <- bp_st
        ds <- res$parameter
        ds["Time"] <- (j-1)*28
        dose_list[[length(dose_list)+1]] <- ds
        #store old info
        old_state <- patient_state$state_id
        # update states
        dbp_obs <- res$bp %>% filter(time==28) #take last observations
        patient_state$bp <- dbp_obs$bp #store observed value
        patient_df$bp0 <- dbp_obs$bp #update system initial conditions
        
        patient_state$discr_bp <- sapply(patient_state$bp, get_discrDBP, 
                                         isInitial=F, simplify = TRUE,
                                         USE.NAMES = TRUE) #discrete value of bp
        
        patient_state$FlagInt[res$parameter$AX==0] <- 1 #patients with AX interr
        patient_state$FlagInt[res$parameter$AX!=0] <- 0 #patients without AX interr
        prev_AX_tmp <- res$parameter$AX
        prev_AH_tmp <- res$parameter$AH 
        prev_AX_tmp[res$parameter$AX==0] <- patient_state$prevAX[res$parameter$AX==0]
        prev_AH_tmp[res$parameter$AX==0] <- patient_state$prevAH[res$parameter$AX==0]
        patient_state$prevAH <- prev_AH_tmp#store previous AH dose
        patient_state$prevAX <- prev_AX_tmp #store previous AX dose
        
        
        patient_state$state_id <- apply(patient_state,1,FUN=get_state_id,
                                        states_map=states_map) #compute new state
        reward <- apply(patient_state,1,FUN=get_reward) #compute reward
        
        
        for (k in 1:Npat) {
            idx_act <- dosing_2_id(patient_state$CurrentAct[k],states_map) # get id action
            
            if (j==24) {
                maxV <- 0
            }else{
                maxV <- max(Q_table_ind.list[[k]][patient_state$state_id[k],]) #maxvalue new state
            }
            
            if (eps<0) {
                # validation step: record reward value
                discRew[k] <- discRew[k] + (reward[k]*(QLparams$gamma^(j-1)))
                
            }else{
                # training step: update Qmatrix
                Q_table_ind.list[[k]][old_state[k],idx_act] <-
                    Q_table_ind.list[[k]][old_state[k],idx_act]+QLparams$alpha*(
                        reward[k]+QLparams$gamma*maxV-Q_table_ind.list[[k]][old_state[k],idx_act]
                    ) 
            }
            
            
            df_k <- patient_state[k,]
            name_action_k <- epsilon_greedy_strategy(df_i=df_k,
                                                     states_action_i = states_action_list[[k]],
                                                     q_i = Q_table_ind.list[[k]],
                                                     eps=eps)
            
            dose_info_k <- parse_dosing(name_action_k)
            
            patient_df$AX[k] <- dose_info_k$AX
            patient_df$AH[k] <- dose_info_k$AH
            if (eps>0) {
                states_action_list[[k]][states_action_list[[k]]$ID==df_k$state_id,name_action_k] <-
                    states_action_list[[k]][states_action_list[[k]]$ID==df_k$state_id,
                                            name_action_k]+1
            }
            patient_state$CurrentAct[k] <- name_action_k
            
        }
        
    }
    
    
    ###### output list
    
    out_list <- list()
    out_list$rew <- discRew
    out_list$bp <- rbindlist(bp_list)
    out_list$dose <- rbindlist(dose_list)
    return(out_list)
    
    
    
    
    
}