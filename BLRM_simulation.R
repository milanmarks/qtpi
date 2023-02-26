# This file implement numerical simulation for BLRM
# Safety rule:
#   1.  if the observed data suggest that there is 25% or higher (posterior) 
#       probability that the DLT rate of a dose is greater than $\delta_2$, 
#       that is, $Pr (p_j > \delta_2 |data) >= 0.25$, 
#       that dose is overdosing and cannot be used to treat patients.
#   2.  Stop the trial if the lowest dose is overdosing

source('code/dosing_scenario.R')
set.seed(654321)
library(rlist)

# Use R package OncoBayes2
library(OncoBayes2)
.user_mc_options <- options(OncoBayes2.MC.warmup=200, OncoBayes2.MC.iter=500, 
                            OncoBayes2.MC.chains=1,
                            OncoBayes2.MC.save_warmup=FALSE)

num_comp <- 1 # one investigational drug
num_inter <- 0 # no drug-drug interactions need to be modeled
num_strata <- 1 # no stratification needed

dref <- 6 # refenrence dose level

path <- "sensitivity_results/weight14/"

# simulation parameters
simN <- 100
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

# weight 8
# W <- rbind(c(0,0.8,1.3,1.8,2.3),
#            c(0,0.5,1,1.5,2),
#            c(0,0,0.25,0.5,1))
# weight matrix for multiple toxicity multiple level
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
  # generate probability of different toxicity level at different dose
  transition_level <- scenario_list[[i]]
  print(transition_level)
  Scn <- Gen_Scn(transition_level = transition_level, 
                 p_r = p_r, p_n = p_n, p_h = p_h, D = 6, W = W, P_T, eps2 = 0.05)
  
  MTD_selection <- numeric()
  patients_allocation <- matrix(NA, nrow=simN, ncol=D)
  p_isotonic <- matrix(NA, nrow=simN, ncol=D)
  FD <- numeric()
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
    
    # clear BLRM model if it exists
    if (exists('blrmfit')) rm(blrmfit)
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
      
      # make decision
      tox_num <- unlist(sapply(DLT, sum))
      record <- tibble(group_id = as.factor('trial_A'),
                            drug_A = 1:D,
                            num_patients = n,
                            num_toxicities= tox_num
      )
      if (!exists('blrmfit')) blrmfit <- blrm_exnex(
        cbind(num_toxicities, num_patients - num_toxicities) ~
          1 + log(drug_A / dref) |
          0 |
          group_id,
        data = record,
        prior_EX_mu_mean_comp = matrix(
          c(logit(1/2),
            log(1)),
          nrow = num_comp,
          ncol = 2
        ),
        prior_EX_mu_sd_comp = matrix(
          c(2,  # sd of intercept
            1), # sd of log-slope
          nrow = num_comp,
          ncol = 2
        ),
        prior_EX_tau_mean_comp = matrix(
          c(0, 0),
          nrow = num_comp,
          ncol = 2
        ),
        prior_EX_tau_sd_comp = matrix(
          c(1, 1),
          nrow = num_comp,
          ncol = 2
        ),
        prior_EX_prob_comp = matrix(1, nrow = num_comp, ncol = 1),
        prior_tau_dist = 0,
        prior_PD = FALSE
      )
      else blrmfit <- update(blrmfit, data = record)
      summ <- summary(blrmfit, interval_prob = c(0, P_T-eps2, P_T+eps2, 1))
      
      # check if first dose is too toxic
      if (summ[1, 8] >= 0.25) {
        mtd <- 0
        break
      }
      
      optional_set <- summ[summ[,8]<0.25,]
      mtd <- as.integer(rownames(optional_set)[which.max(optional_set[,7])])
      if (mtd > d) decision_optimal <- 'E'
      else if (mtd < d) decision_optimal <- 'D'
      else decision_optimal <- 'S'
      
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
      
      # check if failing to de-escalate when 2/3 or >3/6 patients had DLT
      if((sum(DLT[[d_pre]]) >= 2 & length(DLT[[d_pre]]) == 3) | (sum(DLT[[d_pre]]) / length(DLT[[d_pre]]) >= 0.5 & length(DLT[[d_pre]]) >= 6)){
        if (d >= d_pre) {
          fd <- fd + 1
        }
      }
      
      ## stop if the maximum number of patients reached
      if (sum(n) >= samplesize){
        break
      }
      
      ## enroll the next cohort of patients
      cohort_no. <- cohort_no. + 1
    }
    
    # end of one trial
    # determine the MTD
    MTD_selection[sim] <- mtd
    patients_allocation[sim, ] <- n
    FD[sim] <- fd
  }
  
  scenario_abb <- noquote(toupper(paste0(substr(transition_level, 1, 1), collapse = "")))
  write.table(MTD_selection, paste0(path, "BLRM/", scenario_abb, "_MTD_selection.txt"), append=F)
  write.table(patients_allocation, paste0(path, "BLRM/", scenario_abb, "_patients_allocation.txt"), append=F)
  write.table(FD, paste0(path, "BLRM/", scenario_abb, "_FD.txt"), append=F)
}