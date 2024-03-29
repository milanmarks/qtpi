source("code/dosing_scenario.R")
set.seed(654321)
library(rlist)
library(rstan)

################## simulation setup ############################################
path <- 'results/weight8_change_qpt/qtpi1000/'
# simulation parameters
simN <- 1000
qP_T <- 0.2
D <- 6
samplesize <- 27
csize <- 3
startdose <- 1 
eps1 <- .03
eps2 <- .03 
a <- 3
x <- log(1:D/10)
# weight 8
W <- rbind(c(0,0.8,1.3,1.8,2.3),
           c(0,0.5,1,1.5,2),
           c(0,0,0.25,0.5,1))

# weight 12
# W <- rbind(c(0,1,2,4,6),
#            c(0,1,2,3,5),
#            c(0,0,1,2,3))
p_r <- c(0.791, 0.172, 0.032, 0.004, 0.001)
p_n <- c(0.968, 0.029, 0.002, 0.001, 0)
p_h <- c(0.917, 0.070, 0.007, 0.005, 0)

################################################################################
# create equal-length probability intervals 
up.int.bound <- seq(qP_T+eps2, 1, min(eps1+eps2,1-qP_T-eps2)) 
low.int.bound <- sort(-seq(-qP_T+eps1, 0, min(eps1+eps2,qP_T-eps1))) 
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
qLI_index <- which(int.bound < max(low.int.bound))
qEI_index <- which(int.bound >= max(low.int.bound) & int.bound <= min(up.int.bound))[1]
qUI_index <- which(int.bound > min(up.int.bound)) -1
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


for (k in 1:length(scenario_list)){
  transition_level <- scenario_list[[k]]
  print(transition_level)
  Scn <- Gen_Scn(transition_level = transition_level, p_r = p_r, p_n = p_n, p_h = p_h, D = 6, W = W, qP_T, eps2)
  MTD_selection <- numeric()
  patients_allocation <- matrix(NA, nrow=simN, ncol=D)
  EqTP_isotonic <- matrix(NA, nrow=simN, ncol=D)
  FD <- numeric()
  ####################################### Start simulations ###############
  for(sim in 1:simN){
    print(sim)
    n <- rep(0,D)
    d <- startdose
    cohort_no. <- 1
    st <- 0
    nodose <- 0
    toxdose <- D+1
    seldose <- 0
    cohort_sequence <- numeric()
    tox_response <- list(d_1 = list(),
                         d_2 = list(),
                         d_3 = list(),
                         d_4 = list(),
                         d_5 = list(),
                         d_6 = list()) ## list consisting of D=6 doses, subjects treated at each dose, and three toxicity types
    qTP <- list(d_1 = numeric(), 
                d_2 = numeric(), 
                d_3 = numeric(), 
                d_4 = numeric(), 
                d_5 = numeric(), 
                d_6 = numeric()) ## list consisting of D=6 doses, each dose consisting of patients with dTP values
    fd <- 0
    ###############################Single trial ################################ 
    while(st==0){  ## st = 1 indicates the trial must be terminated
      
      cohort_sequence[cohort_no.] <- d  ## track cohort dose allocation
      
      ## generate multi-level toxicity response for three toxicity types
      ## loop over patients in this cohort
      for (i in 1:csize) { 
        tox_i <- matrix(NA, nrow=3, ncol=5)
        p_tox_d <- rbind(Scn$p_all[d,], Scn$p_all[d+D,],  Scn$p_all[d+2*D,])
        tox_i <- t(apply(p_tox_d, 1, rmultinom, n = 1, size = 1))
        tox_response[[d]] <- list.append(tox_response[[d]], tox_i)
      }
      
      ## update qTP at each level
      for(j in 1:max(cohort_sequence)){
        for(i in 1:length(tox_response[[j]])){
          qTP[[j]][i] <- sum(tox_response[[j]][[i]] * W) / sum(apply(W, 1, max))
        }
      }
      
      ## update sample size 
      n[d] <- n[d] + csize   
      
      l <- max(cohort_sequence)
      qtp_sum<- unlist(sapply(qTP, sum))
      data<-list(
        l = l,
        n = as.array(n[1:l]),
        sum_qtp = as.array(qtp_sum[1:l]),
        x = as.array(x[1:l]),
        a = -5,
        b = 5,
        c = 8
      )
      ab_MCMC <- stan('code/qtpi.stan', data = data, refresh = 0)
      params<-extract(ab_MCMC)
      a_post <- params$alpha
      b_post <- params$beta
      EqTP_i <- 1/(1+exp(-a_post-b_post*x[d]))
      Mz_post <- as.numeric(table(cut(EqTP_i, int.bound)))
      Mz_optimal <- which.max(Mz_post)
      decision_optimal <- ifelse(Mz_optimal %in% qLI_index, "E", ifelse(Mz_optimal %in% qEI_index, "S", "D"))
      
      
      ## check safety rule 1 (early termination): if the first dose is too toxic
      EqTP_d1 <- 1/(1+exp(-a_post-b_post*x[1]))
      if(mean(EqTP_d1 > qP_T) > 0.95){
        st <- 1
        break
      }
      
      ## check safety rule 2 (dose exclusion): given the optimal decision is E, is the next dose is too toxic
      if(decision_optimal == "E" & d < toxdose-1){
        EqTP_dnext <- 1/(1+exp(-a_post-b_post*x[d+1]))
        if(mean(EqTP_dnext > qP_T) > 0.95){
          toxdose <- d+1
        }
      }
      
      ## if neither safety rules are revoked, begin dose allocation
      d_pre <- d
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
      if((mean(qTP[[d_pre]]) >= 2/3 & length(qTP[[d_pre]] == 3)) | (mean(qTP[[d_pre]]) >= 0.5) & length(qTP[[d_pre]]) >= 6){
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
    EqTP <- apply(sapply(x, function(x) {1/(1+exp(-a_post-b_post * x))}), 2, mean)
    EqTP_candidate <- EqTP[EqTP <= qP_T + eps2]
    MTD_selection[sim]=ifelse(length(EqTP_candidate)==0,0,which.min(abs(EqTP_candidate - qP_T)))
    patients_allocation[sim, ] <- n
    EqTP_isotonic[sim, ] <- EqTP
    FD[sim] <- fd
  }
  
  ########################## All Trials End, Output Results ####################
  scenario_abb <- noquote(toupper(paste0(substr(transition_level, 1, 1), collapse = "")))
  write.table(MTD_selection, paste0(path, scenario_abb, "_MTD_selection.txt"), append=F)
  write.table(patients_allocation, paste0(path, scenario_abb, "_patients_allocation.txt"), append=F)
  write.table(EqTP_isotonic, paste0(path, scenario_abb, "_EqTP.txt"), append=F)
  write.table(FD, paste0(path, scenario_abb, "_FD.txt"), append=F)
}
