#####------ PREAMBLE ###############
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
cat("\014")
library(dplyr)
library(ggplot2)
library(RsSimulx)
library(survival)
library(data.table)
set.seed(12345)
############## functions ###########################
get_typical_CL <- function(age_gt60,jap,smoke){
  
  typical <- 14.6*(1-0.213*age_gt60)*(1-0.249*jap)*(1+1.02*smoke)
  return(typical)
  
}

get_aucDaily <- function(CLind,fraction,Ddose){
  return(Ddose*fraction/CLind)
  
}

get_logn_samples_BCT <- function(median,cv,lambda,N){
  
  phi <- exp(rnorm(N,mean = 0, sd=cv))
  eta_hat <- ((phi^lambda)-1)/lambda
  
  return(median*exp(eta_hat))
  
  
}


get_ind_CL <- function(age_gt60,jap,smoke,N,sigmaCL=0.599){
  
  typicals <- 14.6*(1-0.213*age_gt60)*(1-0.249*jap)*(1+1.02*smoke)
  ind <- typicals*exp(rnorm(N,mean = 0, sd=sigmaCL))
  return(ind)
  
}


get_bp_ss <- function(base,auc,ahn){
  
  bp_ss <- base*(1+((0.197*auc)/(155.12+auc)))*(1/(1+(0.036*ahn)))
}

############# init list of sampled patients ##############

pat_list <- list()

######### random generation of patients ###################

N<-700000
sigmaAge <- (63-34)/4
u <- runif(N,min =pnorm(-4) ,max =pnorm((80-63)/sigmaAge))
ages <- 63 + sigmaAge*qnorm(u)
min(ages)
max(ages)
#hist(ages)
ages[ages<=60]<-0
ages[ages>60]<-1
jap <- rep(1,N)
smoke <- rep(0,N)
smoke <- runif(N)
smoke[smoke<=0.5] <- 0
smoke[smoke>0.5] <- 1

cl_ind <- get_ind_CL(ages,jap,smoke,N)
dbp0 <- get_logn_samples_BCT(78.9,cv=0.067,lambda=-5.42,N)
dbp0[dbp0>90] <-90 # cirumvent extreme outliers effect
pat_df <- data.frame(id=seq(1,N),CL=cl_ind,BP0=dbp0)

ax_levs <- c(2,3,5,7,10)
ax_names <- c("bid2","bid3","bid5","bid7","bid10")
for (i in 1:length(ax_levs)) {
  a_i <- get_aucDaily(cl_ind,0.457*(1-0.121),ax_levs[i]*2*1e3)
  pat_df[ax_names[i]] <- get_bp_ss(pat_df$BP0,a_i,rep(0,N)) 
}

############## dataframes with 

df_nowindow <- pat_df[pat_df$bid10<90,]

pat_list[[length(pat_list)+1]] <- df_nowindow %>% filter(BP0<70) %>% sample_n(1) %>% 
                                  mutate(Descr="No window <70")
pat_list[[length(pat_list)+1]] <- df_nowindow %>% filter(BP0>=70 & BP0<80) %>% sample_n(1) %>%
                                  mutate(Descr="No window [70-80)")
pat_list[[length(pat_list)+1]] <- df_nowindow %>% filter(BP0>=80 & BP0<=90) %>% sample_n(1) %>%
                                        mutate(Descr="No window [80-90)")

df_hightol <- pat_df[pat_df$bid10>=90 & pat_df$bid10<=100,]

pat_list[[length(pat_list)+1]] <- df_hightol %>% filter(BP0>=70 & BP0<80) %>% sample_n(1) %>%
                                    mutate(Descr="Max Dose optimal [70-80)")
pat_list[[length(pat_list)+1]] <- df_hightol %>% filter(BP0>=80 & BP0<=90) %>% sample_n(1) %>%
                                    mutate(Descr="Max Dose optimal [80-90)")
#############
tox_log <- pat_df
tox_log[,4:8] <-tox_log[,4:8]>=105
tox_log["ToxF"] <- apply(tox_log[,4:8],1,sum)
tox_log_f <- tox_log[tox_log$ToxF>0,]

tname <- list()
for (k in 1:nrow(tox_log_f)) {
  r <- ax_names[which(as.logical(tox_log_f[k,4:8]))]
  tname[[length(tname)+1]] <- r[1]
}
tname <- unlist(tname)
tox_log_f["FirstTox"] <- as.factor(tname)
summary(tox_log_f$FirstTox)
tox_log_f[,4:8] <- pat_df[pat_df$id %in% tox_log_f$id,4:8] 
tox_log_f <- tox_log_f %>% select(!c(ToxF)) %>% 
mutate(Descr=paste("First Dose for severe toxicity:",FirstTox)) 
 
tox_log_f_2 <- tox_log_f %>% filter(FirstTox=="bid2")
pat_list[[length(pat_list)+1]] <- tox_log_f_2 %>% sample_n(1)  %>% select(!c(FirstTox))

tox_log_f_3 <- tox_log_f %>% filter(FirstTox=="bid3")
pat_list[[length(pat_list)+1]] <- tox_log_f_3 %>% sample_n(1)  %>% select(!c(FirstTox))

tox_log_f_5 <- tox_log_f %>% filter(FirstTox=="bid5")
pat_list[[length(pat_list)+1]] <- tox_log_f_5 %>% sample_n(1)  %>% select(!c(FirstTox))

tox_log_f_7 <- tox_log_f %>% filter(FirstTox=="bid7")
pat_list[[length(pat_list)+1]] <- tox_log_f_7 %>% sample_n(1)  %>% select(!c(FirstTox))

tox_log_f_10 <- tox_log_f %>% filter(FirstTox=="bid10")
pat_list[[length(pat_list)+1]] <- tox_log_f_10 %>% sample_n(1)  %>% select(!c(FirstTox))

tox_log_f_10_2ok <- tox_log_f %>% filter(FirstTox=="bid10" & bid2<100)
pat_list[[length(pat_list)+1]] <- tox_log_f_10_2ok %>% sample_n(1)  %>% 
  mutate(Descr="First Dose for severe toxicity: bid10- bid2 target wind")%>%
  select(!c(FirstTox))
##################
no_tox <-  pat_df %>% filter(! id %in% tox_log_f$id)
no_tox[,4:8] <- no_tox[,4:8]>=90 & no_tox[,4:8]<100
no_tox["FlagD"] <- apply(no_tox[,4:8], 1, sum)
no_tox_f <- no_tox[no_tox$FlagD>0,]

tname <- list()
for (k in 1:nrow(no_tox_f)) {
  r <- ax_names[which(as.logical(no_tox_f[k,4:8]))]
  tname[[length(tname)+1]] <- r[length(r)]
}
tname <- unlist(tname)
no_tox_f["FirstDose"] <- as.factor(tname)

summary(no_tox_f$FirstDose)
no_tox_f[,4:8] <- pat_df[pat_df$id %in% no_tox_f$id,4:8] 
no_tox_f <- no_tox_f %>% select(! FlagD) %>%
  filter(FirstDose!="bid10") %>% 
  mutate(Descr=paste("Max window dose:",FirstDose))# those with 10 bid were in df_hightol var.

no_tox_f2 <- no_tox_f %>% filter(FirstDose=="bid2")
pat_list[[length(pat_list)+1]] <- no_tox_f2 %>% sample_n(1) %>% select(! FirstDose)
no_tox_f3 <- no_tox_f %>% filter(FirstDose=="bid3")
pat_list[[length(pat_list)+1]] <- no_tox_f3 %>% sample_n(1) %>% select(! FirstDose)
no_tox_f5 <- no_tox_f %>% filter(FirstDose=="bid5")
pat_list[[length(pat_list)+1]] <- no_tox_f5 %>% sample_n(1) %>% select(! FirstDose)
no_tox_f7 <- no_tox_f %>% filter(FirstDose=="bid7")
pat_list[[length(pat_list)+1]] <- no_tox_f7 %>% sample_n(1) %>% select(! FirstDose)
##################

always_moderate <- pat_df %>% 
  mutate(FlagTox=apply(pat_df[,4:8]>100 & pat_df[,4:8]<105, 1, sum)) %>% filter(FlagTox==5) %>%
  mutate(Descr="Moderate hyper for each Dose")
pat_list[[length(pat_list)+1]] <- always_moderate %>% sample_n(1) %>% select(! FlagTox)

######### final csv ##############

pat_df <- rbindlist(pat_list)
write.csv(pat_df,"patRL.csv", quote = F, row.names = F)