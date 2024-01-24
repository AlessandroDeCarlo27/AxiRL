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
source("test_QL_agents.R")
################### load training results ##########################
load(paste(getwd(),"/training_Qlearning.RData",sep=""))

################## test simulation ###########################################
eps <- -1
test_out <- test_QL_agents(patient_state_initial,patients,eps,states_action_list,
                           Q_table_ind.validation,n_cyc)

################## do some nice plots #########################################
mode <- "png"
if (mode=="pdf") {
    pdf("QLAgents_output.pdf",height=8,width = 9)
    for (i in 1:16) {
        
        
        df<- test_out$bp %>% filter(id==i)
        
        p <- ggplot(data=df) +
            geom_line(aes(x=time,y=bp),size=1.5,linetype="solid")+
            geom_ribbon(aes(x=time,ymin=90,ymax=100),fill="#0072BD",alpha=.2)+
            theme_minimal() + labs(y="Diastolic Blood Pressure [mmHg]",x="Cycle",
                                   title=paste("ID:",i,"- ",patients_raw$Descr[i])) +
            theme(plot.title =  element_text(hjust = 0.5),text = element_text(size=20)) +
            scale_x_continuous(breaks=seq(0,672,28),labels = seq(0,24,1))
        
        
        zz <- test_out$dose %>% filter(id==i) 
        dose_df <- data.frame(dose=c(zz$AX,zz$AH),
                              drug=c(rep("AX",24),rep("AH",24)),
                              time=rep(seq(1,24),2))
        
        
        
        p2 <- ggplot(dose_df, aes(x = time, y = dose))  +
            geom_linerange(
                aes(x = time, ymin = 0, ymax = dose, group = drug), 
                color = "lightgray", size = 1.5,
                position = position_dodge(0.3)
            )+
            geom_point(
                aes(color = drug),
                position = position_dodge(0.3), size = 5
            )+
            scale_color_manual(values = c("#EFC000FF","#0073C2FF"))+
            theme_minimal()+
            theme(plot.title =  element_text(hjust = 0.5),text = element_text(size=20),
                  legend.position="top")+ #+
            labs(y="Administered Amt",x="Cycle")+#+
            scale_y_continuous(breaks=c(0,1,2,3,4,5,7,10),labels=c(0,1,2,3,4,5,7,10)) +#+
            scale_x_continuous(breaks=seq(1,24,1),labels=seq(1,24,1))
        
        grid::grid.draw(gridExtra::grid.arrange(p,p2,ncol = 1, nrow = 2,newpage = TRUE))
    }
    dev.off()
    
    
    
}else{
    
    
    
    for (i in 1:16) {
        
        png(paste0("individualPLots/QLAgents_output",i,".png"))
        df<- test_out$bp %>% filter(id==i)
        
        p <- ggplot(data=df) +
            geom_line(aes(x=time,y=bp),size=1.5,linetype="solid")+
            geom_ribbon(aes(x=time,ymin=90,ymax=100),fill="#0072BD",alpha=.2)+
            theme_minimal() + labs(y="Diastolic Blood Pressure [mmHg]",x="Cycle",
                                   title=paste("ID:",i,"- ",patients_raw$Descr[i])) +
            theme(plot.title =  element_text(hjust = 0.5),text = element_text(size=15)) +
            scale_x_continuous(breaks=seq(0,672,28),labels = seq(0,24,1))
        
        
        zz <- test_out$dose %>% filter(id==i) 
        dose_df <- data.frame(dose=c(zz$AX,zz$AH),
                              drug=c(rep("AX [mg bid]",24),rep("AH [Norm.Dose]",24)),
                              time=rep(seq(1,24),2))
        
        
        
        p2 <- ggplot(dose_df, aes(x = time, y = dose))  +
            geom_linerange(
                aes(x = time, ymin = 0, ymax = dose, group = drug), 
                color = "lightgray", size = 1.5,
                position = position_dodge(0.3)
            )+
            geom_point(
                aes(color = drug),
                position = position_dodge(0.3), size = 5
            )+
            scale_color_manual(values = c("#EFC000FF","#0073C2FF"))+
            theme_minimal()+
            theme(plot.title =  element_text(hjust = 0.5),text = element_text(size=15),
                  legend.position="top")+ #+
            labs(y="Administered Doses",x="Cycle")+#+
            scale_y_continuous(breaks=c(0,1,2,3,4,5,7,10),labels=c(0,1,2,3,4,5,7,10)) +#+
            scale_x_continuous(breaks=seq(1,24,1),labels=seq(1,24,1))
        
        grid::grid.draw(gridExtra::grid.arrange(p,p2,ncol = 1, nrow = 2,newpage = TRUE))
        dev.off()
    }
    
    
    
    
    
}





