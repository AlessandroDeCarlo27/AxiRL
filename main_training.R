########## preamble ##########################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
cat("\014")
library(dplyr)
library(readxl)
library(ggplot2)
library(Matrix)
library(data.table)
library(geomtextpath)
library(RsSimulx)
library(tibble)
source("get_epsilon.R")
source("utils.R")
source("epsilon_greedy_strategy.R")
source("reward.R")
source("train_QL_agents.R")
################## read patient input ###########################

patients_raw <- read.csv("patRL.csv",sep=";")
Npat <- nrow(patients_raw)
patients <- patients_raw %>% mutate(id=1:Npat) %>% 
    select(!c(bid2,bid3,bid5,bid7,bid10,Descr)) %>% rename(bpBase=BP0) 

################ states and actions ##############################

states_map <- read_excel("states.xlsx")
states_action <- read_excel("states_actions_rules.xlsx")
Nstates <- nrow(states_map)
Nactions <- ncol(states_action[,-(1:5)])


############## states action list #################################
states_action_list <- rep(list(states_action),Npat)

############ init Qmatrix #########################################
Q_table_ind.list<-list()
Q_table_ind.list<-lapply(1:Npat,function(x){
    
    Matrix(nrow = Nstates, ncol=Nactions,
           data = 0, sparse = TRUE)
})

Q_table_ind.validation <- Q_table_ind.list

############## scenario params ######################

n_episodes=10#0*1e3 #(training+validation)
n_cyc=24 #number of cycles
freq_validation=2
n_validation_episodes=n_episodes/freq_validation

############# QL params #############################
QLparams <- list()
QLparams$alpha <- 0.1
QLparams$gamma <- 0.95

############### compute initial state #####################
patients$bp0 <- patients$bpBase
patients$AX <- 0
patients$AH <- 0
patient_state_initial <- data.frame(id_pat=patients$id,bp=patients$bp0,
                                    discr_bp = sapply(patients$bp0, get_discrDBP, 
                                                      isInitial=T, simplify = TRUE,
                                                      USE.NAMES = TRUE),
                                    state_id=0,prevAX=0,prevAH=0,FlagInt=0)
patient_state_initial$state_id <- apply(patient_state_initial,1,FUN=get_state_id,
                                        states_map=states_map)
patient_state_initial$CurrentAct <- ""
################### validation buffers ##################
bufMaxRew <- rep(0,Npat)
buf_idx <- rep(0,Npat)
################### TRAINING LOOP #######################


for (episode in 1:n_episodes) {
    
    start.time <- Sys.time()
    #if training
    if (episode%%freq_validation==1) {
        eps<- get_epsilon(ceiling(episode/freq_validation))
        epoch_out <- train_QL_agents(patient_state_initial,patients,eps,states_action_list,
                                     Q_table_ind.list,n_cyc)
        Q_table_ind.list <-epoch_out$Qmat
        states_action_list <- epoch_out$salist
        
    }else{
        
        eps <- -1
        val_out <- train_QL_agents(patient_state_initial,patients,eps,states_action_list,
                                   Q_table_ind.list,n_cyc)
        for (k in 1:Npat) {
            if (val_out$rew[k]>=bufMaxRew[k]) {
                bufMaxRew[k]<-val_out$rew[k]
                Q_table_ind.validation[[k]] <- val_out$Qmat[[k]]
                buf_idx[k] <- episode/2
            }
        }
        
    }
    
    if(episode%%500==0){
        print(paste("episode number: ",episode,sep=""))
        save(list=ls(),file=paste(getwd(),"/training_Qlearning.RData",sep=""))
    }
    
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(paste(as.character(episode),as.character(time.taken),sep=" "))
    
}


save(list=ls(),file=paste(getwd(),"/training_Qlearning.RData",sep=""))

