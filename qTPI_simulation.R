source("code/dosing_scenario.R")
set.seed(654321)
library(rlist)
library(rstan)

################################################################################
## Functions needed
## 1. quasi-likelihood for all patients at all dose levels
total_log_likelihood <- function(a,b,x,qTP){
  y <- sapply(qTP, sum)
  n <- sapply(qTP, length)
  log_lik <- -y*log(1+exp(-a-b*x)) - (n-y)*log(1+exp(a+b*x))
  return(sum(log_lik))
}
## 2. variance for proposal distribution of a
normal_var <- function(x, qTP, a, b){
  y <- sapply(qTP, sum)
  n <- sapply(qTP, length)
  inv_d <- numeric()
  for(d in 1:sum(n!=0)){
    # inv_d[d] <- n[d] / ((1+exp(a+b*x[d]))*(1+exp(-a-b*x[d])))
    inv_d <- n[d]*(y[d]/n[d]+0.005)*(1-y[d]/n[d]+0.005)
  }
  return(1/sum(inv_d))
}
## 3. sample a given b, using MH (prior N(0, 100)); sample b given a using MH (prior unif(0,10))
qTPI_MCMC <- function(qTP, x, iter){
  # empty vector to store samples
  a <- numeric()
  b <- numeric()
  # track acceptance status
  accept_a <- numeric()
  accept_b <- numeric()
  # initial value for a, b
  a[1] <- runif(n=1, min=-5, max=5)
  b[1] <- runif(n=1, min=0, max=8)
  accept_a[1] <- 1
  accept_b[1] <- 1
  # initial proposal distribution
  a_current <- a[1]
  b_current <- b[1]
  # MCMC begins
  for(t in 2:iter){
    # 1. given current value of b, propose new a
    # a_propose <- rnorm(n=1, mean=a_current, sd=sqrt(normal_var(x=x, qTP=qTP, a=a_current, b=b_current)))
    a_propose <- runif(n=1, min=-5, max=5)
    # calculate acceptance status for proposed a
    log_ratio_a <- total_log_likelihood(a=a_propose, b=b_current,x=x,qTP=qTP) - total_log_likelihood(a=a_current, b=b_current,x=x,qTP=qTP)  
    
    # accept a?
    if (log_ratio_a > log(runif(1))){
      a_current <- a_propose
      accept_a[t] <- 1
    }else{
      a_current <- a_current
      accept_a[t] <- 0
    }
    a[t] <- a_current
    
    # 2. given current value of a, propose new b
    b_propose <- runif(n=1, min=0, max=8)
    # calculate acceptance ratio
    log_ratio_b <- total_log_likelihood(a=a_current, b=b_propose,x=x,qTP=qTP) - total_log_likelihood(a=a_current, b=b_current,x=x,qTP=qTP)
    # accept b?
    if (log_ratio_b > log(runif(1))){
      b_current <- b_propose
      accept_b[t] <- 1
    }else{
      b_current <- b_current
      accept_b[t] <- 0
    }
    b[t] <- b_current
  }
  return(list(b_posterior = b[(0.3*iter):iter],
              a_posterior = a[(0.3*iter):iter],
              accept_a = accept_a,
              accept_b = accept_b))
}

## 4. pava is the pool adjacent violator algorithm to perform isotonic transformation for the posterior means later
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
## 5. betavar computes variances of beta distributions 
betavar<-function(a,b){
  resp <- a*b/((a+b)^2*(a+b+1))
  return(resp)
}

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
