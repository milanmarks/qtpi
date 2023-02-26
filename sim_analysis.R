library(ggplot2)
library(kableExtra)
library(dplyr)
source("code/dosing_scenario.R")

path = "results/weight8/"

# weight 8
W <- rbind(c(0,0.8,1.3,1.8,2.3),
           c(0,0.5,1,1.5,2),
           c(0,0,0.25,0.5,1))
# weight 12
# W <- rbind(c(0,1,2,4,6),
#            c(0,1,2,3,5),
#            c(0,0,1,2,3))

P_T <- 0.3
eps2 <- 0.05
p_r <- c(0.791, 0.172, 0.032, 0.004, 0.001)
p_n <- c(0.968, 0.029, 0.002, 0.001, 0)
p_h <- c(0.917, 0.070, 0.007, 0.005, 0)
MTD_choice <- "MTD_qTP"
record <- record_qTP(p_r,p_n,p_h,W,P_T,eps2)
MTD_by_qTP <- record$MTD_by_qTP
MTD_by_p <- record$MTD_by_p
# write.csv(format(record$qTP, digits=2), file=paste0(path, 'qTP.csv'))
# write.csv(format(record$p, digits = 2), file=paste0(path, 'p.csv'))

scn_list <- c("LLL", "LLM", "LLH", "LML", "LMM", "LMH", "LHL", "LHM", "LHH",
                  "MLL", "MLM", "MLH", "MML", "MMM", "MMH", "MHL", "MHM", "MHH",
                  "HLL", "HLM", "HLH", "HML", "HMM", "HMH", "HHL", "HHM", "HHH")

################################################################################
## safety
################################################################################
## percent selection above MTD

safety_summary <- data.frame(scenario = scn_list,
                             MTD_qTP = MTD_by_qTP,
                             MTD_p = MTD_by_p,
                             POS_qTPI = rep(NA, 27),
                             POS_mTPI2 = rep(NA, 27),
                             POS_BOIN = rep(NA, 27),
                             POS_CRM = rep(NA, 27),
                             POS_BLRM = rep(NA, 27),
                             POA_qTPI = rep(NA, 27),
                             POA_mTPI2 = rep(NA, 27),
                             POA_BOIN = rep(NA, 27),
                             POA_CRM = rep(NA, 27),
                             POA_BLRM = rep(NA, 27))
# qTPI
for(i in 1:length(scn_list)){
  scn <- scn_list[i]
  filename_MTD <- paste0(path, "qtpi/", scn, "_MTD_selection.txt")
  filename_patients <- paste0(path, "qtpi/", scn, "_patients_allocation.txt")
  MTD_selection <- read.table(filename_MTD, header = T)
  patients_allocation <- read.table(filename_patients, header = T)
  patients_allocation$POA <- NULL
  safety_summary[i, "scenario"] <- scn
  safety_summary[i, "POS_qTPI"] <- mean(MTD_selection > safety_summary[i, MTD_choice])
  for(n in 1:nrow(patients_allocation)){
    if(safety_summary[i, MTD_choice] < 6){
      patients_allocation[n, "POA"] <- sum(patients_allocation[n, (safety_summary[i, MTD_choice]+1):6]) / sum(patients_allocation[n,1:6])
    }else{
      patients_allocation[n, "POA"] <- 0
    }
  }
  safety_summary[i, "POA_qTPI"] <- mean(patients_allocation$POA)
}

# mTPI2
for(i in 1:length(scn_list)){
  scn <- scn_list[i]
  filename_MTD <- paste0(path, "mTPI2/", scn, "_MTD_selection.txt")
  filename_patients <- paste0(path, "mTPI2/", scn, "_patients_allocation.txt")
  MTD_selection <- read.table(filename_MTD, header = T)
  patients_allocation <- read.table(filename_patients, header = T)
  patients_allocation$POA <- NULL
  safety_summary[i, "scenario"] <- scn
  safety_summary[i, "POS_mTPI2"] <- mean(MTD_selection > safety_summary[i, MTD_choice])
  for(n in 1:nrow(patients_allocation)){
    if(safety_summary[i, MTD_choice] < 6){
      patients_allocation[n, "POA"] <- sum(patients_allocation[n, (safety_summary[i, MTD_choice]+1):6]) / sum(patients_allocation[n,1:6])
    }else{
      patients_allocation[n, "POA"] <- 0
    }
  }
  safety_summary[i, "POA_mTPI2"] <- mean(patients_allocation$POA)
}

# BOIN
for(i in 1:length(scn_list)){
  scn <- scn_list[i]
  filename_MTD <- paste0(path, "BOIN/", scn, "_MTD_selection.txt")
  filename_patients <- paste0(path, "BOIN/", scn, "_patients_allocation.txt")
  MTD_selection <- read.table(filename_MTD, header = T)
  patients_allocation <- read.table(filename_patients, header = T)
  patients_allocation$POA <- NULL
  safety_summary[i, "scenario"] <- scn
  safety_summary[i, "POS_BOIN"] <- mean(MTD_selection > safety_summary[i, MTD_choice])
  for(n in 1:nrow(patients_allocation)){
    if(safety_summary[i, MTD_choice] < 6){
      patients_allocation[n, "POA"] <- sum(patients_allocation[n, (safety_summary[i, MTD_choice]+1):6]) / sum(patients_allocation[n,1:6])
    }else{
      patients_allocation[n, "POA"] <- 0
    }
  }
  safety_summary[i, "POA_BOIN"] <- mean(patients_allocation$POA)
}

# CRM
for(i in 1:length(scn_list)){
  scn <- scn_list[i]
  filename_MTD <- paste0(path, "CRM/", scn, "_MTD_selection.txt")
  filename_patients <- paste0(path, "CRM/", scn, "_patients_allocation.txt")
  MTD_selection <- read.table(filename_MTD, header = T)
  patients_allocation <- read.table(filename_patients, header = T)
  patients_allocation$POA <- NULL
  safety_summary[i, "scenario"] <- scn
  safety_summary[i, "POS_CRM"] <- mean(MTD_selection > safety_summary[i, MTD_choice])
  for(n in 1:nrow(patients_allocation)){
    if(safety_summary[i, MTD_choice] < 6){
      patients_allocation[n, "POA"] <- sum(patients_allocation[n, (safety_summary[i, MTD_choice]+1):6]) / sum(patients_allocation[n,1:6])
    }else{
      patients_allocation[n, "POA"] <- 0
    }
  }
  safety_summary[i, "POA_CRM"] <- mean(patients_allocation$POA)
}

# BLRM
for(i in 1:length(scn_list)){
  scn <- scn_list[i]
  filename_MTD <- paste0(path, "BLRM/", scn, "_MTD_selection.txt")
  filename_patients <- paste0(path, "BLRM/", scn, "_patients_allocation.txt")
  MTD_selection <- read.table(filename_MTD, header = T)
  patients_allocation <- read.table(filename_patients, header = T)
  patients_allocation$POA <- NULL
  safety_summary[i, "scenario"] <- scn
  safety_summary[i, "POS_BLRM"] <- mean(MTD_selection > safety_summary[i, MTD_choice])
  for(n in 1:nrow(patients_allocation)){
    if(safety_summary[i, MTD_choice] < 6){
      patients_allocation[n, "POA"] <- sum(patients_allocation[n, (safety_summary[i, MTD_choice]+1):6]) / sum(patients_allocation[n,1:6])
    }else{
      patients_allocation[n, "POA"] <- 0
    }
  }
  safety_summary[i, "POA_BLRM"] <- mean(patients_allocation$POA)
}

## Automatically generate latex code for supplemental table
safety_summary[, c(1,7,8,6,5,4,12,13,11,10,9)] %>%
  kable(
    caption = "Full Simulation Results Comparing the Safety Performance 
    of CRM, BLRM, BOIN, mTPI2 and qTPI designs under 27 Dosing Scenarios",
    format = "latex",
    escape = F,
    digits = c(3,3,3,3,3,3,3,3,3,3),
    row.names = NA,
    col.names = c("Scenario", "POS: CRM", "POS: BLRM", "POS: BOIN", 
                  "POS: mTPI2", "POS: qTPI", "POA: CRM", "POA: BLRM",
                  "POA: BOIN", "POA: mTPI2", "POA: qTPI"),
    align = "lccccccccccc")%>%
  kable_styling(bootstrap_options = "condensed", font_size = 7) %>%
  column_spec(1, bold = T) %>%
  collapse_rows(columns = 1, valign = "top")%>%
  footnote(number = c("Scenario ``LMH'' represents: low, moderate and high renal, neurological and hematological toxicities"))


## visualization
safety_summary$Scn <- 1:27
safety_summary <- safety_summary[order(safety_summary$POS_qTPI), ]
safety_summary$id <- 1:27
safety_summary_long1 <- reshape(
  safety_summary[, c("Scn", "id", "POS_qTPI", "POS_mTPI2", "POS_BOIN", "POS_CRM", "POS_BLRM")],
  direction = "long",
  varying = c("POS_qTPI", "POS_mTPI2", "POS_BOIN", "POS_CRM", "POS_BLRM"),
  v.names = "value",
  idvar = c("Scn", "id"),
  timevar = "Method",
  times = c("qTPI", "mTPI2", "BOIN", "CRM", "BLRM"))
safety_summary_long1$Metric <- "POS"
row.names(safety_summary_long1) <- NULL

ggplot(safety_summary_long1, aes(x=id, y=value, linetype=Method, color=Method))+
  geom_line(aes(linetype=Method), lwd=0.8)+
  scale_x_continuous(breaks=1:27,labels=as.character(safety_summary$Scn))+
  theme_minimal()+
  scale_linetype_manual(name="",values=c('dashed','dashed','dashed','dashed','solid')) +
  scale_color_manual(name="", values=c('red','blue','forestgreen','orange','black')) +
  labs(x="Dosing Scenario", y="Percent of Overdose Selection", title="")+
  theme(axis.text.x = element_text(colour = 'black', size = 12,hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(colour = 'black', size = 12),
        axis.title.x =element_text(size=15),
        axis.title.y =element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15, ),
        legend.key.size = unit(1, "cm"),
        legend.position="top") 

ggsave("pos.png", path = paste0(path, "graphs/"), 
       bg = "white", width = 8, height = 5)

safety_summary <- safety_summary[order(safety_summary$POA_qTPI),]
safety_summary$id <- 1:27
safety_summary_long2 <- reshape(
  safety_summary[, c("Scn", "id", "POA_qTPI", "POA_mTPI2", "POA_BOIN", 
                     "POA_CRM", "POA_BLRM")],
  direction = "long",
  varying = c("POA_qTPI", "POA_mTPI2", "POA_BOIN", "POA_CRM", "POA_BLRM"),
  v.names = "value",
  idvar = c("Scn", "id"),
  timevar = "Method",
  times = c("qTPI", "mTPI2", "BOIN", "CRM", "BLRM"))
safety_summary_long2$Metric <- "POA"
row.names(safety_summary_long2) <- NULL

ggplot(safety_summary_long2, aes(x=id, y=value, linetype=Method))+
  geom_line(aes(color=Method), lwd=0.8)+
  scale_linetype_manual(name="",values=c('dashed','dashed','dashed','dashed','solid')) +
  scale_color_manual(name="", values=c('red','blue','forestgreen','orange','black')) +
  scale_x_continuous(breaks=1:27,labels=as.character(safety_summary$Scn))+
  theme_minimal()+
  labs(x="Dosing Scenario", y="Percent of Overdose Allocation", title="")+
  theme(axis.text.x = element_text(colour = 'black', size = 12,hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(colour = 'black', size = 12),
        axis.title.x =element_text(size=15),
        axis.title.y =element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.key.size = unit(1, "cm"),
        legend.position="top") 

ggsave("poa.png", path = paste0(path, "graphs/"), 
       bg = "white", width = 8, height = 5)

safety_summary_MTD <- data.frame(Method = rep(NA, 10),
                                 MTD6 = rep(NA, 10),
                                 MTD5 = rep(NA, 10),
                                 MTD4 = rep(NA, 10),
                                 MTD3 = rep(NA, 10))
MTD <- c(6, 5, 4, 3)
safety_summary_MTD$Method <- c("POS: qTPI", "POS: mTPI2", "POS: BOIN", "POS: CRM", "POS: BLRM",
                               "POA: qTPI", "POA: mTPI2", "POA: BOIN", "POS: CRM", "POS: BLRM")
for(i in 1:4){
  safety_summary_MTD[1, i+1] <- mean(safety_summary[safety_summary$MTD_qTP == MTD[i], "POS_qTPI"])
  safety_summary_MTD[2, i+1] <- mean(safety_summary[safety_summary$MTD_qTP == MTD[i], "POS_mTPI2"])
  safety_summary_MTD[3, i+1] <- mean(safety_summary[safety_summary$MTD_qTP == MTD[i], "POS_BOIN"])
  safety_summary_MTD[4, i+1] <- mean(safety_summary[safety_summary$MTD_qTP == MTD[i], "POS_CRM"])
  safety_summary_MTD[5, i+1] <- mean(safety_summary[safety_summary$MTD_qTP == MTD[i], "POS_BLRM"])
  safety_summary_MTD[6, i+1] <- mean(safety_summary[safety_summary$MTD_qTP == MTD[i], "POA_qTPI"])
  safety_summary_MTD[7, i+1] <- mean(safety_summary[safety_summary$MTD_qTP == MTD[i], "POA_mTPI2"])
  safety_summary_MTD[8, i+1] <- mean(safety_summary[safety_summary$MTD_qTP == MTD[i], "POA_BOIN"])
  safety_summary_MTD[9, i+1] <- mean(safety_summary[safety_summary$MTD_qTP == MTD[i], "POA_CRM"])
  safety_summary_MTD[10, i+1] <- mean(safety_summary[safety_summary$MTD_qTP == MTD[i], "POA_BLRM"])
}

################################################################################
## accuracy
################################################################################
## percent selection above MTD
accuracy_summary <- data.frame(scenario = scn_list,
                               MTD_qTP = MTD_by_qTP,
                               MTD_p = MTD_by_p,
                               PCS_qTPI = rep(NA, 27),
                               PCS_mTPI2 = rep(NA, 27),
                               PCS_BOIN = rep(NA, 27),
                               PCS_CRM = rep(NA, 27),
                               PCS_BLRM = rep(NA, 27),
                               PCA_qTPI = rep(NA, 27),
                               PCA_mTPI2 = rep(NA, 27),
                               PCA_BOIN = rep(NA, 27),
                               PCA_CRM = rep(NA, 27),
                               PCA_BLRM = rep(NA, 27))
# qTPI
for(i in 1:length(scn_list)){
  scn <- scn_list[i]
  filename_MTD <- paste0(path, "qtpi/", scn, "_MTD_selection.txt")
  filename_patients <- paste0(path, "qtpi/", scn, "_patients_allocation.txt")
  MTD_selection <- read.table(filename_MTD, header = T)
  patients_allocation <- read.table(filename_patients, header = T)
  patients_allocation$PCA <- NULL
  accuracy_summary[i, "scenario"] <- scn
  accuracy_summary[i, "PCS_qTPI"] <- mean(MTD_selection == accuracy_summary[i,MTD_choice])
  for(n in 1:nrow(patients_allocation)){
    patients_allocation[n, "PCA"] <- patients_allocation[n, accuracy_summary[i,MTD_choice]] / sum(patients_allocation[n,1:6])
  }
  accuracy_summary[i, "PCA_qTPI"] <- mean(patients_allocation$PCA)
}

# mTPI2
for(i in 1:length(scn_list)){
  scn <- scn_list[i]
  filename_MTD <- paste0(path, "mtpi2/", scn, "_MTD_selection.txt")
  filename_patients <- paste0(path, "mtpi2/", scn, "_patients_allocation.txt")
  MTD_selection <- read.table(filename_MTD, header = T)
  patients_allocation <- read.table(filename_patients, header = T)
  patients_allocation$PCA <- NULL
  accuracy_summary[i, "scenario"] <- scn
  accuracy_summary[i, "PCS_mTPI2"] <- mean(MTD_selection == accuracy_summary[i,MTD_choice])
  for(n in 1:nrow(patients_allocation)){
    patients_allocation[n, "PCA"] <- patients_allocation[n, accuracy_summary[i,MTD_choice]] / sum(patients_allocation[n,1:6])
  }
  accuracy_summary[i, "PCA_mTPI2"] <- mean(patients_allocation$PCA)
}

# BOIN
for(i in 1:length(scn_list)){
  scn <- scn_list[i]
  filename_MTD <- paste0(path, "BOIN/", scn, "_MTD_selection.txt")
  filename_patients <- paste0(path, "BOIN/", scn, "_patients_allocation.txt")
  MTD_selection <- read.table(filename_MTD, header = T)
  patients_allocation <- read.table(filename_patients, header = T)
  patients_allocation$PCA <- NULL
  accuracy_summary[i, "scenario"] <- scn
  accuracy_summary[i, "PCS_BOIN"] <- mean(MTD_selection == accuracy_summary[i,MTD_choice])
  for(n in 1:nrow(patients_allocation)){
    patients_allocation[n, "PCA"] <- patients_allocation[n, accuracy_summary[i,MTD_choice]] / sum(patients_allocation[n,1:6])
  }
  accuracy_summary[i, "PCA_BOIN"] <- mean(patients_allocation$PCA)
}

# CRM
for(i in 1:length(scn_list)){
  scn <- scn_list[i]
  filename_MTD <- paste0(path, "CRM/", scn, "_MTD_selection.txt")
  filename_patients <- paste0(path, "CRM/", scn, "_patients_allocation.txt")
  MTD_selection <- read.table(filename_MTD, header = T)
  patients_allocation <- read.table(filename_patients, header = T)
  patients_allocation$PCA <- NULL
  accuracy_summary[i, "scenario"] <- scn
  accuracy_summary[i, "PCS_CRM"] <- mean(MTD_selection == accuracy_summary[i,MTD_choice])
  for(n in 1:nrow(patients_allocation)){
    patients_allocation[n, "PCA"] <- patients_allocation[n, accuracy_summary[i,MTD_choice]] / sum(patients_allocation[n,1:6])
  }
  accuracy_summary[i, "PCA_CRM"] <- mean(patients_allocation$PCA)
}

# BLRM
for(i in 1:length(scn_list)){
  scn <- scn_list[i]
  filename_MTD <- paste0(path, "BLRM/", scn, "_MTD_selection.txt")
  filename_patients <- paste0(path, "BLRM/", scn, "_patients_allocation.txt")
  MTD_selection <- read.table(filename_MTD, header = T)
  patients_allocation <- read.table(filename_patients, header = T)
  patients_allocation$PCA <- NULL
  accuracy_summary[i, "scenario"] <- scn
  accuracy_summary[i, "PCS_BLRM"] <- mean(MTD_selection == accuracy_summary[i,MTD_choice])
  for(n in 1:nrow(patients_allocation)){
    patients_allocation[n, "PCA"] <- patients_allocation[n, accuracy_summary[i,MTD_choice]] / sum(patients_allocation[n,1:6])
  }
  accuracy_summary[i, "PCA_BLRM"] <- mean(patients_allocation$PCA)
}

## Automatically generate latex code for supplemental table
accuracy_summary[, c(1,7,8,6,5,4,12,13,11,10,9)] %>%
  kable(
    caption = "Full Simulation Results Comparing the Accuracy Performance of 
    CRM, BLRM, BOIN, mTPI2 and qTPI designs under 27 Dosing Scenarios",
    format = "latex",
    escape = F,
    digits = c(3,3,3,3,3,3,3,3,3,3),
    row.names = NA,
    col.names = c("Scenario", "PCS: CRM", "PCS: BLRM", "PCS: BOIN", "PCS: mTPI2",
                  "PCS: qTPI", "PCA: CRM", "PCA: BLRM", "PCA: BOIN", 
                  "PCA: mTPI2", "PCA: qTPI"),
    align = "lcccccccccc")%>%
  kable_styling(bootstrap_options = "condensed", font_size = 7) %>%
  column_spec(1, bold = T) %>%
  collapse_rows(columns = 1, valign = "top")%>%
  footnote(number = c("Scenario ``LMH'' represents: low, moderate and high renal, neurological and hematological toxicities"))

## visualization
accuracy_summary$Scn <- 1:27
accuracy_summary <- accuracy_summary[order(accuracy_summary$PCS_qTPI), ]
accuracy_summary$id <- 1:27
accuracy_summary_long1 <- reshape(
  accuracy_summary[, c("Scn", "id", "PCS_qTPI", "PCS_mTPI2", "PCS_BOIN", 
                       "PCS_CRM", "PCS_BLRM")],
  direction = "long",
  varying = c("PCS_qTPI", "PCS_mTPI2", "PCS_BOIN", "PCS_CRM", "PCS_BLRM"),
  v.names = "value",
  idvar = c("Scn", "id"),
  timevar = "Method",
  times = c("qTPI", "mTPI2", "BOIN", "CRM", "BLRM"))
accuracy_summary_long1$Metric <- "PCS"
row.names(accuracy_summary_long1) <- NULL

ggplot(accuracy_summary_long1, aes(x=id, y=value, linetype=Method, color=Method))+
  geom_line(aes(linetype=Method), lwd=0.8)+
  scale_linetype_manual(name="",values=c('dashed','dashed','dashed','dashed','solid')) +
  scale_color_manual(name="", values=c('red','blue','forestgreen','orange','black')) +
  scale_x_continuous(breaks=1:27,labels=as.character(accuracy_summary$Scn))+
  theme_minimal()+
  labs(x="Dosing Scenario", y="Percent of Correct Selection", title="", linetype="Design")+
  theme(axis.text.x = element_text(colour = 'black', size = 12,hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(colour = 'black', size = 12),
        axis.title.x =element_text(size=15),
        axis.title.y =element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.key.size = unit(1, "cm"),
        legend.position="top") 
ggsave("pcs.png", path = paste0(path, "graphs/"), 
       bg = "white", width = 8, height = 5)

accuracy_summary <- accuracy_summary[order(accuracy_summary$PCA_qTPI), ]
accuracy_summary$id <- 1:27
accuracy_summary_long2 <- reshape(
  accuracy_summary[, c("Scn", "id", "PCA_qTPI", "PCA_mTPI2", "PCA_BOIN", 
                       "PCA_CRM", "PCA_BLRM")],
  direction = "long",
  varying = c("PCA_qTPI", "PCA_mTPI2", "PCA_BOIN", "PCA_CRM", "PCA_BLRM"),
  v.names = "value",
  idvar = c("Scn", "id"),
  timevar = "Method",
  times = c("qTPI", "mTPI2", "BOIN", "CRM", "BLRM"))
accuracy_summary_long2$Metric <- "PCA"
row.names(accuracy_summary_long2) <- NULL

ggplot(accuracy_summary_long2, aes(x=id, y=value))+
  geom_line(aes(linetype=Method, color=Method), lwd=0.8)+
  scale_linetype_manual(name="",values=c('dashed','dashed','dashed','dashed','solid')) +
  scale_color_manual(name="", values=c('red','blue','forestgreen','orange','black')) +
  scale_x_continuous(breaks=1:27,labels=as.character(accuracy_summary$Scn))+
  theme_minimal()+
  labs(x="Dosing Scenario", y="Percent of Correct Allocation", title="", linetype="Design")+
  theme(axis.text.x = element_text(colour = 'black', size = 12,hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(colour = 'black', size = 12),
        axis.title.x =element_text(size=15),
        axis.title.y =element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.key.size = unit(1, "cm"),
        legend.position="top") 
ggsave("pca.png", path = paste0(path, "graphs/"), 
       bg = "white", width = 8, height = 5)

accuracy_summary_MTD <- data.frame(Method = rep(NA, 10),
                                   MTD6 = rep(NA, 10),
                                   MTD5 = rep(NA, 10),
                                   MTD4 = rep(NA, 10),
                                   MTD3 = rep(NA, 10))
MTD <- c(6, 5, 4, 3)
accuracy_summary_MTD$Method <- c("PCS: qTPI", "PCS: mTPI2", "PCS: BOIN", "PCS: CRM", "PCS: BLRM",
                                 "PCA: qTPI", "PCA: mTPI2", "PCA: BOIN", "PCA: CRM", "PCA: BLRM")
for(i in 1:4){
  accuracy_summary_MTD[1, i+1] <- mean(accuracy_summary[accuracy_summary$MTD_qTP == MTD[i], "PCS_qTPI"])
  accuracy_summary_MTD[2, i+1] <- mean(accuracy_summary[accuracy_summary$MTD_qTP == MTD[i], "PCS_mTPI2"])
  accuracy_summary_MTD[3, i+1] <- mean(accuracy_summary[accuracy_summary$MTD_qTP == MTD[i], "PCS_BOIN"])
  accuracy_summary_MTD[4, i+1] <- mean(accuracy_summary[accuracy_summary$MTD_qTP == MTD[i], "PCS_CRM"])
  accuracy_summary_MTD[5, i+1] <- mean(accuracy_summary[accuracy_summary$MTD_qTP == MTD[i], "PCS_BLRM"])
  accuracy_summary_MTD[6, i+1] <- mean(accuracy_summary[accuracy_summary$MTD_qTP == MTD[i], "PCA_qTPI"])
  accuracy_summary_MTD[7, i+1] <- mean(accuracy_summary[accuracy_summary$MTD_qTP == MTD[i], "PCA_mTPI2"])
  accuracy_summary_MTD[8, i+1] <- mean(accuracy_summary[accuracy_summary$MTD_qTP == MTD[i], "PCA_BOIN"])
  accuracy_summary_MTD[9, i+1] <- mean(accuracy_summary[accuracy_summary$MTD_qTP == MTD[i], "PCA_CRM"])
  accuracy_summary_MTD[10, i+1] <- mean(accuracy_summary[accuracy_summary$MTD_qTP == MTD[i], "PCA_BLRM"])
}

################################################################################
## reliability
################################################################################
## overdosing more than 50% of the patients
reliability_summary <- data.frame(scenario = scn_list,
                                  MTD_qTP = MTD_by_qTP,
                                  MTD_p = MTD_by_p,
                                  OP50_qTPI = rep(NA, 27),
                                  OP50_mTPI2 = rep(NA, 27),
                                  OP50_BOIN = rep(NA, 27),
                                  OP50_CRM = rep(NA, 27),
                                  OP50_BLRM = rep(NA, 27),
                                  FD_qTPI = rep(NA, 27),
                                  FD_mTPI2 = rep(NA, 27),
                                  FD_BOIN = rep(NA, 27),
                                  FD_CRM = rep(NA, 27),
                                  FD_BLRM = rep(NA, 27))
# qTPI
for(i in 1:length(scn_list)){
  scn <- scn_list[i]
  filename_FD <- paste0(path, "qtpi/", scn, "_FD.txt")
  filename_patients <- paste0(path, "qtpi/", scn, "_patients_allocation.txt")
  patients_allocation <- read.table(filename_patients, header = T)
  failD <- read.table(filename_FD, header = T)
  patients_allocation$OP50 <- NULL
  reliability_summary[i, "scenario"] <- scn
  for(n in 1:nrow(patients_allocation)){
    if(reliability_summary[i, "MTD_qTP"] < 6){
      patients_allocation[n, "OP50"] <- ifelse(sum(patients_allocation[n, (reliability_summary[i,"MTD_qTP"]+1):6]) >= 0.5*sum(patients_allocation[n,1:6]), 1, 0)
    }else{
      patients_allocation[n, "OP50"] <- 0
    }
  }
  reliability_summary[i, "OP50_qTPI"] <- mean(patients_allocation$OP50)
  reliability_summary[i, "FD_qTPI"] <- mean(as.integer(as.logical(failD[,1])))
}


# mTPI2
for(i in 1:length(scn_list)){
  scn <- scn_list[i]
  filename_FD <- paste0(path, "mtpi2/", scn, "_FD.txt")
  filename_patients <- paste0(path, "mtpi2/", scn, "_patients_allocation.txt")
  patients_allocation <- read.table(filename_patients, header = T)
  failD <- read.table(filename_FD, header = T)
  patients_allocation$OP50 <- NULL
  reliability_summary[i, "scenario"] <- scn
  for(n in 1:nrow(patients_allocation)){
    if(reliability_summary[i, "MTD_qTP"] < 6){
      patients_allocation[n, "OP50"] <- ifelse(sum(patients_allocation[n, (reliability_summary[i,"MTD_qTP"]+1):6]) >= 0.5*sum(patients_allocation[n,1:6]), 1, 0)
    }else{
      patients_allocation[n, "OP50"] <- 0
    }
  }
  reliability_summary[i, "OP50_mTPI2"] <- mean(patients_allocation$OP50)
  reliability_summary[i, "FD_mTPI2"] <- mean(as.integer(as.logical(failD[,1])))
}

# BOIN
for(i in 1:length(scn_list)){
  scn <- scn_list[i]
  filename_FD <- paste0(path, "BOIN/", scn, "_FD.txt")
  filename_patients <- paste0(path, "BOIN/", scn, "_patients_allocation.txt")
  patients_allocation <- read.table(filename_patients, header = T)
  failD <- read.table(filename_FD, header = T)
  patients_allocation$OP50 <- NULL
  reliability_summary[i, "scenario"] <- scn
  for(n in 1:nrow(patients_allocation)){
    if(reliability_summary[i, "MTD_qTP"] < 6){
      patients_allocation[n, "OP50"] <- ifelse(sum(patients_allocation[n, (reliability_summary[i,"MTD_qTP"]+1):6]) >= 0.5*sum(patients_allocation[n,1:6]), 1, 0)
    }else{
      patients_allocation[n, "OP50"] <- 0
    }
  }
  reliability_summary[i, "OP50_BOIN"] <- mean(patients_allocation$OP50)
  reliability_summary[i, "FD_BOIN"] <- mean(as.integer(as.logical(failD[,1])))
}

# CRM 
for(i in 1:length(scn_list)){
  scn <- scn_list[i]
  filename_FD <- paste0(path, "CRM/", scn, "_FD.txt")
  filename_patients <- paste0(path, "CRM/", scn, "_patients_allocation.txt")
  patients_allocation <- read.table(filename_patients, header = T)
  failD <- read.table(filename_FD, header = T)
  patients_allocation$OP50 <- NULL
  reliability_summary[i, "scenario"] <- scn
  for(n in 1:nrow(patients_allocation)){
    if(reliability_summary[i, "MTD_qTP"] < 6){
      patients_allocation[n, "OP50"] <- ifelse(sum(patients_allocation[n, (reliability_summary[i,"MTD_qTP"]+1):6]) >= 0.5*sum(patients_allocation[n,1:6]), 1, 0)
    }else{
      patients_allocation[n, "OP50"] <- 0
    }
  }
  reliability_summary[i, "OP50_CRM"] <- mean(patients_allocation$OP50)
  reliability_summary[i, "FD_CRM"] <- mean(as.integer(as.logical(failD[,1])))
}

# BLRM
for(i in 1:length(scn_list)){
  scn <- scn_list[i]
  filename_FD <- paste0(path, "BLRM/", scn, "_FD.txt")
  filename_patients <- paste0(path, "BLRM/", scn, "_patients_allocation.txt")
  patients_allocation <- read.table(filename_patients, header = T)
  failD <- read.table(filename_FD, header = T)
  patients_allocation$OP50 <- NULL
  reliability_summary[i, "scenario"] <- scn
  for(n in 1:nrow(patients_allocation)){
    if(reliability_summary[i, "MTD_qTP"] < 6){
      patients_allocation[n, "OP50"] <- ifelse(sum(patients_allocation[n, (reliability_summary[i,"MTD_qTP"]+1):6]) >= 0.5*sum(patients_allocation[n,1:6]), 1, 0)
    }else{
      patients_allocation[n, "OP50"] <- 0
    }
  }
  reliability_summary[i, "OP50_BLRM"] <- mean(patients_allocation$OP50)
  reliability_summary[i, "FD_BLRM"] <- mean(as.integer(as.logical(failD[,1])))
}

reliability_summary[, c(1,7,8,6,5,4,12,13,11,10,9)] %>%
  kable(
    caption = "Full Simulation Results Comparing the Reliability Performance 
    of CRM, BLRM, BOIN, mTPI2 and qTPI designs under 27 Dosing Scenarios",
    format = "latex",
    escape = F,
    digits = c(3,3,3,3,3,3,3,3,3,3),
    row.names = NA,
    col.names = c("Scenario", "OP50: CRM", "OP50: BLRM", "OP50: BOIN", 
                  "OP50: mTPI2","OP50: qTPI", "FD: CRM", "FD: BLRM",
                  "FD: BOIN", "FD: mTPI2", "FD: qTPI"),
    align = "lcccccccccc")%>%
  kable_styling(bootstrap_options = "condensed", font_size = 7) %>%
  column_spec(1, bold = T) %>%
  collapse_rows(columns = 1, valign = "top")%>%
  footnote(number = c("Scenario ``LMH'' represents: low, moderate and high renal, neurological and hematological toxicities"))


## visualization
reliability_summary$Scn <- 1:27
reliability_summary <- reliability_summary[order(reliability_summary$OP50_qTPI), ]
reliability_summary$id <- 1:27
reliability_summary_long <- reshape(
reliability_summary[, c("Scn", "id", "OP50_qTPI", "OP50_mTPI2", "OP50_BOIN", 
                          "OP50_CRM", "OP50_BLRM")],
  direction = "long",
  varying = c("OP50_qTPI", "OP50_mTPI2", "OP50_BOIN", "OP50_CRM", "OP50_BLRM"),
  v.names = "value",
  idvar = c("Scn", "id"),
  timevar = "Method",
  times = c("qTPI", "mTPI2", "BOIN", "CRM", "BLRM"))
reliability_summary_long$Metric <- "OP50"
row.names(reliability_summary_long) <- NULL

ggplot(reliability_summary_long, aes(x=id, y=value))+
  geom_line(aes(linetype=Method, color=Method), lwd=0.8)+
  scale_linetype_manual(name="",values=c('dashed','dashed','dashed','dashed','solid')) +
  scale_color_manual(name="", values=c('red','blue','forestgreen','orange','black')) +
  scale_x_continuous(breaks=1:27,labels=as.character(reliability_summary$Scn))+
  theme_minimal()+
  labs(x="Dosing Scenario", y="Percent of Overdosing 50% Patients", title="", linetype="Design")+
  theme(axis.text.x = element_text(colour = 'black', size = 12,hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(colour = 'black', size = 12),
        axis.title.x =element_text(size=15),
        axis.title.y =element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.key.size = unit(1, "cm"),
        legend.position="top") 
ggsave("op50.png", path = paste0(path, "graphs/"), 
       bg = "white", width = 8, height = 5)

reliability_summary$Scn <- 1:27
reliability_summary <- reliability_summary[order(reliability_summary$FD_qTPI), ]
reliability_summary$id <- 1:27
reliability_summary_long2 <- reshape(
reliability_summary[, c("Scn", "id", "FD_qTPI", "FD_mTPI2", "FD_BOIN", 
                          "FD_CRM", "FD_BLRM")],
  direction = "long",
  varying = c("FD_qTPI", "FD_mTPI2", "FD_BOIN", "FD_CRM", "FD_BLRM"),
  v.names = "value",
  idvar = c("Scn", "id"),
  timevar = "Method",
  times = c("qTPI", "mTPI2", "BOIN", "CRM", "BLRM"))
reliability_summary_long$Metric <- "FD"
row.names(reliability_summary_long) <- NULL

ggplot(reliability_summary_long2, aes(x=id, y=value))+
  geom_line(aes(linetype=Method, color=Method), lwd=0.8)+
  scale_linetype_manual(name="",values=c('dashed','dashed','dashed','dashed','solid')) +
  scale_color_manual(name="", values=c('red','blue','forestgreen','orange','black')) +
  scale_x_continuous(breaks=1:27,labels=as.character(reliability_summary$Scn))+
  theme_minimal()+
  labs(x="Dosing Scenario", y="Percentage of fail to deescalate", title="", linetype="Design")+
  theme(axis.text.x = element_text(colour = 'black', size = 12,hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(colour = 'black', size = 12),
        axis.title.x =element_text(size=15),
        axis.title.y =element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.key.size = unit(1, "cm"),
        legend.position="top") 
ggsave("FD.png", path = paste0(path, "graphs/"), 
       bg = "white", width = 8, height = 5)

reliability_summary_MTD <- data.frame(Method = rep(NA, 10),
                                      MTD6 = rep(NA, 10),
                                      MTD5 = rep(NA, 10),
                                      MTD4 = rep(NA, 10),
                                      MTD3 = rep(NA, 10))
MTD <- c(6, 5, 4, 3)
reliability_summary_MTD$Method <- c("OP50: qTPI", "OP50: mTPI2", "OP50: BOIN", "OP50: CRM", "OP50: BLRM",
                                    "FD: qTPI", "FD: mTPI2", "FD: BOIN", "FD: CRM", "FD: BLRM")
for(i in 1:4){
  reliability_summary_MTD[1, i+1] <- mean(reliability_summary[reliability_summary$MTD_qTP == MTD[i], "OP50_qTPI"])
  reliability_summary_MTD[2, i+1] <- mean(reliability_summary[reliability_summary$MTD_qTP == MTD[i], "OP50_mTPI2"])
  reliability_summary_MTD[3, i+1] <- mean(reliability_summary[reliability_summary$MTD_qTP == MTD[i], "OP50_BOIN"])
  reliability_summary_MTD[4, i+1] <- mean(reliability_summary[reliability_summary$MTD_qTP == MTD[i], "OP50_CRM"])
  reliability_summary_MTD[5, i+1] <- mean(reliability_summary[reliability_summary$MTD_qTP == MTD[i], "OP50_BLRM"])
  reliability_summary_MTD[6, i+1] <- mean(reliability_summary[reliability_summary$MTD_qTP == MTD[i], "FD_qTPI"])
  reliability_summary_MTD[7, i+1] <- mean(reliability_summary[reliability_summary$MTD_qTP == MTD[i], "FD_mTPI2"])
  reliability_summary_MTD[8, i+1] <- mean(reliability_summary[reliability_summary$MTD_qTP == MTD[i], "FD_BOIN"])
  reliability_summary_MTD[9, i+1] <- mean(reliability_summary[reliability_summary$MTD_qTP == MTD[i], "FD_CRM"])
  reliability_summary_MTD[10, i+1] <- mean(reliability_summary[reliability_summary$MTD_qTP == MTD[i], "FD_BLRM"])
}

################################################################################
## early stop and average sample size
################################################################################
## percent selection above MTD
other_summary <- data.frame(scenario = scn_list,
                            MTD_qTP = MTD_by_qTP,
                            MTD_p = MTD_by_p,
                            SampleSize_qTPI = rep(NA, 27),
                            SampleSize_mTPI2 = rep(NA, 27),
                            SampleSize_BOIN = rep(NA, 27),
                            SampleSize_CRM = rep(NA, 27),
                            SampleSize_BLRM = rep(NA, 27),
                            EarlyStop_qTPI = rep(NA, 27),
                            EarlyStop_mTPI2 = rep(NA, 27),
                            EarlyStop_BOIN = rep(NA, 27),
                            EarlyStop_CRM = rep(NA, 27),
                            EarlyStop_BLRM = rep(NA, 27))
# qTPI
for(i in 1:length(scn_list)){
  scn <- scn_list[i]
  filename_patients <- paste0(path, "qtpi/", scn, "_patients_allocation.txt")
  patients_allocation <- read.table(filename_patients, header = T)
  patients_allocation$EarlyStop <- NULL
  other_summary[i, "scenario"] <- scn
  other_summary[i, "SampleSize_qTPI"] <- mean(rowSums(patients_allocation))
  for(n in 1:nrow(patients_allocation)){
    patients_allocation[n, "EarlyStop"] <- sum(patients_allocation[n, 1:6]) < 27
  }
  other_summary[i, "EarlyStop_qTPI"] <- mean(patients_allocation$EarlyStop)
}


# mTPI2
for(i in 1:length(scn_list)){
  scn <- scn_list[i]
  filename_patients <- paste0(path, "mtpi2/", scn, "_patients_allocation.txt")
  patients_allocation <- read.table(filename_patients, header = T)
  patients_allocation$EarlyStop <- NULL
  other_summary[i, "scenario"] <- scn
  other_summary[i, "SampleSize_mTPI2"] <- mean(rowSums(patients_allocation))
  for(n in 1:nrow(patients_allocation)){
    patients_allocation[n, "EarlyStop"] <- sum(patients_allocation[n, 1:6]) < 27
  }
  other_summary[i, "EarlyStop_mTPI2"] <- mean(patients_allocation$EarlyStop)
}

# BOIN
for(i in 1:length(scn_list)){
  scn <- scn_list[i]
  filename_patients <- paste0(path, "BOIN/", scn, "_patients_allocation.txt")
  patients_allocation <- read.table(filename_patients, header = T)
  patients_allocation$EarlyStop <- NULL
  other_summary[i, "scenario"] <- scn
  other_summary[i, "SampleSize_BOIN"] <- mean(rowSums(patients_allocation))
  for(n in 1:nrow(patients_allocation)){
    patients_allocation[n, "EarlyStop"] <- sum(patients_allocation[n, 1:6]) < 27
  }
  other_summary[i, "EarlyStop_BOIN"] <- mean(patients_allocation$EarlyStop)
}

# CRM
for(i in 1:length(scn_list)){
  scn <- scn_list[i]
  filename_patients <- paste0(path, "CRM/", scn, "_patients_allocation.txt")
  patients_allocation <- read.table(filename_patients, header = T)
  patients_allocation$EarlyStop <- NULL
  other_summary[i, "scenario"] <- scn
  other_summary[i, "SampleSize_CRM"] <- mean(rowSums(patients_allocation))
  for(n in 1:nrow(patients_allocation)){
    patients_allocation[n, "EarlyStop"] <- sum(patients_allocation[n, 1:6]) < 27
  }
  other_summary[i, "EarlyStop_CRM"] <- mean(patients_allocation$EarlyStop)
}

# BLRM
for(i in 1:length(scn_list)){
  scn <- scn_list[i]
  filename_patients <- paste0(path, "BLRM/", scn, "_patients_allocation.txt")
  patients_allocation <- read.table(filename_patients, header = T)
  patients_allocation$EarlyStop <- NULL
  other_summary[i, "scenario"] <- scn
  other_summary[i, "SampleSize_BLRM"] <- mean(rowSums(patients_allocation))
  for(n in 1:nrow(patients_allocation)){
    patients_allocation[n, "EarlyStop"] <- sum(patients_allocation[n, 1:6]) < 27
  }
  other_summary[i, "EarlyStop_BLRM"] <- mean(patients_allocation$EarlyStop)
}

other_summary[, c(1,7,8,6,5,4,12,13,11,10,9)] %>%
  kable(
    caption = "Full Simulation Results Comparing the Sample Size and Early Stop 
    of CRM, BLRM, BOIN, mTPI2 and qTPI designs under 27 Dosing Scenarios",
    format = "latex",
    escape = F,
    digits = c(3,3,3,3,3,3,3,3,3,3),
    row.names = NA,
    col.names = c("Scenario", "SZ: CRM","SZ: BLRM","SZ: BOIN","SZ: mTPI2", 
                  "SZ: qTPI", "ES: CRM", "ES: BLRM", "ES: BOIN", "ES: mTPI2", 
                  "ES: qTPI"),
    align = "lccccccccccc")%>%
  kable_styling(bootstrap_options = "condensed", font_size = 7) %>%
  column_spec(1, bold = T) %>%
  collapse_rows(columns = 1, valign = "top")%>%
  footnote(number = c("Scenario ``LMH'' represents: low, moderate and high renal, neurological and hematological toxicities"))


# visualization
other_summary$Scn <- 1:27
other_summary <- other_summary[order(other_summary$EarlyStop_mTPI2), ]
other_summary$id <- 1:27
other_summary_long1 <- reshape(
  other_summary[, c("Scn", "id", "SampleSize_qTPI", "SampleSize_mTPI2", 
                    "SampleSize_BOIN", "SampleSize_CRM", "SampleSize_BLRM")],
  direction = "long",
  varying = c("SampleSize_qTPI", "SampleSize_mTPI2", "SampleSize_BOIN", 
              "SampleSize_CRM", "SampleSize_BLRM"),
  v.names = "value",
  idvar = c("Scn", "id"),
  timevar = "Method",
  times = c("qTPI", "mTPI2", "BOIN", "CRM", "BLRM"))
other_summary_long1$Metric <- "Sample Size"
row.names(other_summary_long1) <- NULL

other_summary_long2 <- reshape(
  other_summary[, c("Scn", "id", "EarlyStop_qTPI", "EarlyStop_mTPI2", 
                    "EarlyStop_BOIN", "EarlyStop_CRM", "EarlyStop_BLRM")],
  direction = "long",
  varying = c("EarlyStop_qTPI", "EarlyStop_mTPI2", "EarlyStop_BOIN", 
              "EarlyStop_CRM", "EarlyStop_BLRM"),
  v.names = "value",
  idvar = c("Scn", "id"),
  timevar = "Method",
  times = c("qTPI", "mTPI2", "BOIN", "CRM", "BLRM"))
other_summary_long2$Metric <- "Early Stop"
row.names(other_summary_long2) <- NULL

other_summary_long <- rbind(other_summary_long1, other_summary_long2)

ggplot(other_summary_long2, aes(x=id, y=value))+
  geom_line(aes(linetype=Method, color=Method), lwd=0.8)+
  scale_linetype_manual(name="",values=c('dashed','dashed','dashed','dashed','solid')) +
  scale_color_manual(name="", values=c('red','blue','forestgreen','orange','black')) +
  scale_x_continuous(breaks=1:27,labels=as.character(other_summary$Scn))+
  theme_minimal()+
  labs(x="Dosing Scenario", y="Percent of Early Stop", title="", linetype="Design")+
  theme(axis.text.x = element_text(colour = 'black', size = 12,hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(colour = 'black', size = 12),
        axis.title.x =element_text(size=15),
        axis.title.y =element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.key.size = unit(1, "cm"),
        legend.position="top") 
ggsave("earlystop.png", path = paste0(path, "graphs/"), 
       bg = "white", width = 8, height = 5)

other_summary_MTD <- data.frame(Method = rep(NA, 10),
                                MTD6 = rep(NA, 10),
                                MTD5 = rep(NA, 10),
                                MTD4 = rep(NA, 10),
                                MTD3 = rep(NA, 10))
MTD <- c(6, 5, 4, 3)
other_summary_MTD$Method <- c("SampleSize: qTPI", "SampleSize: mTPI2", 
                              "SampleSize: BOIN", "SampleSize: CRM", 
                              "SampleSize: BLRM",
                              "EarlyStop: qTPI", "EarlyStop: mTPI2", 
                              "EarlyStop: BOIN", "EarlyStop: CRM",
                              "EarlyStop: BLRM")
for(i in 1:4){
  other_summary_MTD[1, i+1] <- mean(other_summary[other_summary$MTD_qTP == MTD[i], "SampleSize_qTPI"])
  other_summary_MTD[2, i+1] <- mean(other_summary[other_summary$MTD_qTP == MTD[i], "SampleSize_mTPI2"])
  other_summary_MTD[3, i+1] <- mean(other_summary[other_summary$MTD_qTP == MTD[i], "SampleSize_BOIN"])
  other_summary_MTD[4, i+1] <- mean(other_summary[other_summary$MTD_qTP == MTD[i], "SampleSize_CRM"])
  other_summary_MTD[5, i+1] <- mean(other_summary[other_summary$MTD_qTP == MTD[i], "SampleSize_BLRM"])
  other_summary_MTD[6, i+1] <- mean(other_summary[other_summary$MTD_qTP == MTD[i], "EarlyStop_qTPI"])
  other_summary_MTD[7, i+1] <- mean(other_summary[other_summary$MTD_qTP == MTD[i], "EarlyStop_mTPI2"])
  other_summary_MTD[8, i+1] <- mean(other_summary[other_summary$MTD_qTP == MTD[i], "EarlyStop_BOIN"])
  other_summary_MTD[9, i+1] <- mean(other_summary[other_summary$MTD_qTP == MTD[i], "EarlyStop_CRM"])
  other_summary_MTD[10, i+1] <- mean(other_summary[other_summary$MTD_qTP == MTD[i], "EarlyStop_BLRM"])
}

all_summary_MTD <- format(
  rbind(safety_summary_MTD, accuracy_summary_MTD, 
        reliability_summary_MTD, other_summary_MTD),
  digits = 4,
  nsmall = 4
)
write.csv(all_summary_MTD, file = 'results/weight8/all_sum.csv')



