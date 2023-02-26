# This file implement numerical simulation for CRM
# Two rules are imposed
# 1.  No skipping dose level
# 2.  Stop the trial if the posterior probability that the DLT probability of 
#     the lowest dose is greater than target exceeds a certain threshold, 
#     such as 0.95, that is, P(p1 > target | data) > 0.95

source('code/dosing_scenario.R')
set.seed(654321)
library(rlist)
library(dfcrm)

path <- 'sensitivity_results/weight14/'

# simulation parameters
simN <- 1000  
# target DLT rate
P_T <- 0.3
eps2 <- 0.05
# Number of Doses
D <- 6
# overall samplesize
samplesize <- 27
# cohort size
csize <- 3
startdose <- 1 

# weight matrix for multiple toxicity multiple level
# W <- rbind(c(0,0.8,1.3,1.8,2.3),
#            c(0,0.5,1,1.5,2),
#            c(0,0,0.25,0.5,1))
W <- rbind(c(0,2,3,4,5),
           c(0,1,2,3,4),
           c(0,0,0.5,1,2))

# initial probability vector of grade 0 to grade 4 renal, 
# neurological and hematological toxicities at dose level 1
p_r <- c(0.791, 0.172, 0.032, 0.004, 0.001)
p_n <- c(0.968, 0.029, 0.002, 0.001, 0)
p_h <- c(0.917, 0.070, 0.007, 0.005, 0)

# low, moderate, high transition probability matrix for dose level increase
trans_matrix_L <- as.matrix(rbind(c(0.8,0.2,0,0,0),
                                  c(0,0.8,0.2,0,0),
                                  c(0,0,0.8,0.2,0),
                                  c(0,0,0,0.8,0.2),
                                  c(0,0,0,0,1)))
trans_matrix_M <- as.matrix(rbind(c(0.6,0.4,0,0,0),
                                  c(0,0.6,0.4,0,0),
                                  c(0,0,0.6,0.4,0),
                                  c(0,0,0,0.6,0.4),
                                  c(0,0,0,0,1)))
trans_matrix_H <- as.matrix(rbind(c(0.4,0.6,0,0,0),
                                  c(0,0.4,0.6,0,0),
                                  c(0,0,0.4,0.6,0),
                                  c(0,0,0,0.4,0.6),
                                  c(0,0,0,0,1)))
# list of different scenario types
scenario_list <- list(c("low","low","low"),
                      c("low","low","moderate"),
                      c("low","low","high"),
                      c("low","moderate","low"),
                      c("low","moderate","moderate"),
                      c("low","moderate","high"),
                      c("low","high","low"),
                      c("low","high","moderate"),
                      c("low","high","high"),
                      c("moderate","low","low"),
                      c("moderate","low","moderate"),
                      c("moderate","low","high"),
                      c("moderate","moderate","low"),
                      c("moderate","moderate","moderate"),
                      c("moderate","moderate","high"),
                      c("moderate","high","low"),
                      c("moderate","high","moderate"),
                      c("moderate","high","high"),
                      c("high","low","low"),
                      c("high","low","moderate"),
                      c("high","low","high"),
                      c("high","moderate","low"),
                      c("high","moderate","moderate"),
                      c("high","moderate","high"),
                      c("high","high","low"),
                      c("high","high","moderate"),
                      c("high","high","high"))


for(i in 1:length(scenario_list)){
  # generate scenario
  transition_level <- scenario_list[[i]]
  Scn <- Gen_Scn(transition_level = transition_level, 
                 p_r = p_r, p_n = p_n, p_h = p_h, D = 6, W = W, P_T, eps2 = 0.05)
  
  MTD_selection <- numeric()
  patients_allocation <- matrix(NA, nrow=simN, ncol=D)
  p_isotonic <- matrix(NA, nrow=simN, ncol=D)
  FD <- numeric()
  
  # get prior based on current scenario for CRM model
  prior_MTD <- Scn$MTD_by_qTP
  prior <- getprior(halfwidth = 0.06, target = P_T, nu = prior_MTD, nlevel = D)
  
  for(sim in 1:simN){
    print(sim)
    n <- rep(0,D)
    d <- startdose
    cohort_no. <- 1
    st <- 0
    toxdose <- D+1
    cohort_sequence <- numeric()
    
    # list consisting of D=6 doses, 
    # subjects treated at each dose, and three toxicity types
    tox_response <- list(d_1 = list(),
                         d_2 = list(),
                         d_3 = list(),
                         d_4 = list(),
                         d_5 = list(),
                         d_6 = list()) 
    # list consisting of D=6 doses, 
    # each dose consisting of patients with binary DLT status
    DLT <- list(d_1 = numeric(),
                d_2 = numeric(),
                d_3 = numeric(),
                d_4 = numeric(),
                d_5 = numeric(),
                d_6 = numeric())
    fd <- 0
    
    while(st==0){  ## st = 1 indicates the trial must be terminated
      
      cohort_sequence[cohort_no.] <- d  ## track cohort dose allocation
      
      ## generate multi-level toxicity response for three toxicity types
      for(i in 1:csize){ ## loop over patients in this cohort
        tox_i <- matrix(NA, nrow=3, ncol=5)
        for(j in 1:3){ ## loop over toxicity type
          tox_i[j, ] <- rmultinom(n=1, size=1, prob = Scn$p_all[d+(j-1)*D,])
        }
        tox_response[[d]] <- list.append(tox_response[[d]], tox_i)
      }
      
      ## update binary DLT at each level
      for(j in 1:max(cohort_sequence)){
        for(i in 1:length(tox_response[[j]])){
          DLT[[j]][i] <- as.numeric((tox_response[[j]][[i]][1,4] == 1) | (tox_response[[j]][[i]][1,5] ==1) | (tox_response[[j]][[i]][2,4] ==1) | (tox_response[[j]][[i]][2,5] ==1) | (tox_response[[j]][[i]][3,5] ==1))
        }
      }
      
      ## update sample size 
      n[d] <- n[d] + csize   
      
      # make decision
      tox <- unlist(DLT)
      level <- rep(seq(n), n)
      crm_model <- crm(prior, target = P_T, tox, level)
      mtd <- crm_model$mtd
      if (d < mtd){
        decision_optimal <- 'E'
      } else if (d == mtd){
        decision_optimal <- 'S'
      } else{
        decision_optimal <- 'D'
      }
      
      est <- crm_model$estimate
      var <- crm_model$post.var
      
      ## check safety rule
      if(d == 1){
        if( pnorm(log(P_T / prior[1]), est, var, lower.tail = F) > 0.95 ){
          st <- 1
          break
        }
      }
      
      # check if failing to de-escalate when 2/3 or >3/6 patients had DLT
      if((sum(DLT[[d]]) >= 2 & length(DLT[[d]]) == 3) | (sum(DLT[[d]]) / length(DLT[[d]]) >= 0.5 & length(DLT[[d]]) >= 6)){
        fd <- fd + ifelse(decision_optimal == "D", 0, 1)
      }
      
      ## if neither safety rules are revoked, begin dose allocation
      if(decision_optimal == "D"){
        if(d != 1){
          d <- d-1
        }else{
          d <- 1
        }
      }else if(decision_optimal == "E"){
        if(d != D){
          d <- min(d+1, toxdose-1)
        }else{
          d <- D
        }
      }else{
        d <- min(d, toxdose-1)
      }
      ## stop if the maximum number of patients reached
      if (sum(n) >= samplesize){
        break
      }
      ## enroll the next cohort of patients
      cohort_no. <- cohort_no. + 1
    } # end of one cohort
    
    ## determine the MTD from the single trial
    MTD_selection[sim] <- mtd
    patients_allocation[sim, ] <- n
    FD[sim] <- fd
  }
  
  scenario_abb <- noquote(toupper(paste0(substr(transition_level, 1, 1), collapse = "")))
  write.table(MTD_selection, paste0(path, "CRM/", scenario_abb, "_MTD_selection.txt"), append=F)
  write.table(patients_allocation, paste0(path, "CRM/", scenario_abb, "_patients_allocation.txt"), append=F)
  write.table(FD, paste0(path, "CRM/", scenario_abb, "_FD.txt"), append=F)
}
