source("code/dosing_scenario.R")
set.seed(654321)
library(rlist)

################################################################################
## Functions needed
## 1. pava is the pool adjacent violator algorithm to perform isotonic transformation for the posterior means later
pava <- function (x, wt = rep(1, length(x))) 
{
  n <- length(x)
  if (n <= 1) 
    return(x)
  if (any(is.na(x)) || any(is.na(wt))) {
    stop("Missing values in 'x' or 'wt' not allowed")
  }
  lvlsets <- (1:n)
  repeat {
    viol <- (as.vector(diff(x)) < 0)
    if (!(any(viol))) 
      break
    i <- min((1:(n - 1))[viol])
    lvl1 <- lvlsets[i]
    lvl2 <- lvlsets[i + 1]
    ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
    x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
    lvlsets[ilvl] <- lvl1
  }
  x
}
## 2. betavar computes variances of beta distributions 
betavar<-function(a,b){
  resp <- a*b/((a+b)^2*(a+b+1))
  return(resp)
}

################## simulation setup ############################################
path <- 'sensitivity_results/weight14/mtpi2/'
# simulation parameters
simN <- 1000  
P_T <- 0.3
D <- 6
samplesize <- 27
csize <- 3
startdose <- 1 
eps1 <- .05
eps2 <- .05 
a <- 1
b <- 1
W <- rbind(c(0,2,3,4,5),
           c(0,1,2,3,4),
           c(0,0,0.5,1,2))
p_r <- c(0.791, 0.172, 0.032, 0.004, 0.001)
p_n <- c(0.968, 0.029, 0.002, 0.001, 0)
p_h <- c(0.917, 0.070, 0.007, 0.005, 0)

################################################################################
# create equal-length probability intervals 
up.int.bound <- seq(P_T+eps2, 1, min(eps1+eps2,1-P_T-eps2)) 
low.int.bound <- sort(-seq(-P_T+eps1, 0, min(eps1+eps2,P_T-eps1))) 
if (low.int.bound[1] != 0){
  if (up.int.bound[length(up.int.bound)] != 1)
    int.bound = c(0,low.int.bound,up.int.bound,1)
  else 
    int.bound = c(0,low.int.bound,up.int.bound)
} else {
  if (up.int.bound[length(up.int.bound)] != 1)
    int.bound = c(low.int.bound,up.int.bound,1)
  else 
    int.bound = c(low.int.bound,up.int.bound)
}
LI_index <- which(int.bound < max(low.int.bound))
EI_index <- which(int.bound >= max(low.int.bound) & int.bound <= min(up.int.bound))[1]
UI_index <- which(int.bound > min(up.int.bound)) -1

################################################################################
# set up scenarios
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

################################################################################
## run simulation
################################################################################
for(i in 1:length(scenario_list)){  
  transition_level <- scenario_list[[i]]
  Scn <- Gen_Scn(transition_level = transition_level, p_r = p_r, p_n = p_n, p_h = p_h, D = 6, W = W, P_T, eps2)
  MTD_selection <- numeric()
  patients_allocation <- matrix(NA, nrow=simN, ncol=D)
  p_isotonic <- matrix(NA, nrow=simN, ncol=D)
  FD <- numeric()
  ####################################### Start simulations ###############
  for(sim in 1:simN){
    print(sim)
    n <- rep(0,D)
    d <- startdose
    cohort_no. <- 1
    st <- 0
    toxdose <- D+1
    cohort_sequence <- numeric()
    tox_response <- list(d_1 = list(),
                         d_2 = list(),
                         d_3 = list(),
                         d_4 = list(),
                         d_5 = list(),
                         d_6 = list()) ## list consisting of D=6 doses, subjects treated at each dose, and three toxicity types
    
    DLT <- list(d_1 = numeric(),
                d_2 = numeric(),
                d_3 = numeric(),
                d_4 = numeric(),
                d_5 = numeric(),
                d_6 = numeric()) ## list consisting of D=6 doses, each dose consisting of patients with binary DLT status
    fd <- 0
    ###############################Single trial ################################ 
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
      
      ## update the posterior of Mz
      Mz_post <- numeric()
      for(z in 1:(length(int.bound)-1)){
        Mz_post[z] <- pbeta(
          int.bound[z+1], 1 + sum(DLT[[d]]),1 + length(DLT[[d]])-sum(DLT[[d]])
          )-pbeta(
            int.bound[z], 1 + sum(DLT[[d]]), 1 + length(DLT[[d]])-sum(DLT[[d]])
            )
        Mz_post[z] <- Mz_post[z] / (int.bound[z+1] - int.bound[z])
      }
      Mz_optimal <- which.max(Mz_post)
      decision_optimal <- ifelse(Mz_optimal %in% LI_index, "E", ifelse(Mz_optimal %in% EI_index, "S", "D"))
      
      ## check safety rule 1 (early termination): if the first dose is too toxic
      if(d == 1){
        if(1-pbeta(P_T, 0.005+sum(DLT[[1]]), 0.005+length(DLT[[1]])-sum(DLT[[1]])) > 0.95){
          st <- 1
          break
        }
      }
      
      ## check safety rule 2 (dose exclusion): given the optimal decision is E, is the next dose is too toxic
      ## then exclude the next and higher dose and stay at the current dose
      if(decision_optimal == "E" & d < toxdose-1){
        if(1-pbeta(P_T, 0.005+sum(DLT[[d+1]]), 0.005+length(DLT[[d+1]])-sum(DLT[[d+1]])) > 0.95){
          toxdose <- d+1
        }
      }
      d_pre <- d
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
      
      ## check if failing to de-escalate when 2/3 or >3/6 patients had DLT
      if((sum(DLT[[d_pre]]) >= 2 & length(DLT[[d_pre]]) == 3) | (sum(DLT[[d_pre]]) / length(DLT[[d_pre]]) >= 0.5 & length(DLT[[d_pre]]) >= 6)){
        fd <- fd + ifelse(decision_optimal == "D", 0, 1)
      }
      
      ## enroll the next cohort of patients
      cohort_no. <- cohort_no. + 1
      
      ## stop if the maximum number of patients reached
      if (sum(n) >= samplesize){
        st <- 1
        break
      }
    }
    
    ############################### Single Trial Summary #######################
    ## determine the MTD from the single trial
    ## first, update posterior of b using the full data observed which is the results from the last iteration
    ## second, calculate the posterior mean and variance of EqTP_d at each dose level
    p_mean <- numeric()
    p_var <- numeric()
    for(d in 1:D){
      p_mean[d] <- (sum(DLT[[d]])+0.005) / (length(DLT[[d]])+0.01)
      p_var[d] <- betavar(sum(DLT[[d]])+0.005, length(DLT[[d]])-sum(DLT[[d]])+0.005)
    }
    ## third, apply istonic regression to the posterior mean of EqTP_d
    p_mean_isotonic <- pava(x=p_mean, wt=1/p_var)
    ## break tie
    for(d in 1:D){
      p_mean_isotonic[d] <- p_mean_isotonic[d] + (d-1)*1E-10
    }
    ## fourth, select the MTD as the highest dose with min|EqTP - qP_T| and EqTP <= 0.35
    p_candidate <- p_mean_isotonic[p_mean_isotonic <= P_T + eps2]
    if (length(p_candidate)){
      MTD_selection[sim] <- which.min(abs(p_candidate - P_T))  
    } else{
      MTD_selection[sim] <- 1
    }
    
    patients_allocation[sim, ] <- n
    p_isotonic[sim, ] <- p_mean_isotonic
    FD[sim] <- fd
  }
  
  ########################## All Trials End, Output Results ####################
  scenario_abb <- noquote(toupper(paste0(substr(transition_level, 1, 1), collapse = "")))
  write.table(MTD_selection, paste0(path, scenario_abb, "_MTD_selection.txt"), append=F)
  write.table(patients_allocation, paste0(path, scenario_abb, "_patients_allocation.txt"), append=F)
  write.table(p_isotonic, paste0(path, scenario_abb, "_p.txt"), append=F)
  write.table(FD, paste0(path, scenario_abb, "_FD.txt"), append=F)
}

