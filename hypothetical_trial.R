################################################################################
## Hypothetical Trial Example
## qTPI - logistic design
################################################################################
rm(list = ls())
set.seed(654321)
library(rlist)
library(rstan)


# define mtpi2 function
mtpi2 <- function(d, DLT, tox_response) {
  
  # update the DLT matrix 
  for (i in (length(tox_response[[d]])-2):length(tox_response[[d]])){
    flag <- 0
    for (j in 1:6){
      if (tox_response[[d]][[i]][j, 4] == 1 | tox_response[[d]][[i]][j, 5] == 1)
        flag <- 1
    }
    DLT[[d]] <- list.append(DLT[[d]], flag)
  }
  
  # prior distribution of dose i toxicity probability $p_i$
  # $p_i \sim Beta(alpha, beta)$, default to be $Beta(1, 1)$ 
  alpha <- 1
  beta <- 1
  
  # posterior distribution of $p_i$
  # $p_i\mid data \sim Beta(alpha + n_d, beta + n_d - x_d)$
  post_alpha <- alpha + sum(DLT[[d]])
  post_beta <- beta + length(DLT[[d]]) - sum(DLT[[d]])
  
  Mz_post_mtpi2 <- numeric()
  for(z in 1:(length(int.bound)-1)){
    Mz_post_mtpi2[z] <- pbeta(
      int.bound[z+1], 
      post_alpha,
      post_beta
    )-
      pbeta(
        int.bound[z], 
        post_alpha,
        post_beta
      )
  }
  Mz_optimal_mtpi2 <- which.max(Mz_post_mtpi2)
  decision_optimal_mtpi2 <- ifelse(
    Mz_optimal_mtpi2 %in% qLI_index, 
    "E", 
    ifelse(Mz_optimal_mtpi2 %in% qEI_index, "S", "D")
    )
  return (list(
    DLT = DLT,
    Mz_post_mtpi2 = Mz_post_mtpi2,
    Mz_optimal_mtpi2 = Mz_optimal_mtpi2,
    decision_optimal_mtpi2 = decision_optimal_mtpi2
  ))
}

x <- log(c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
qP_T <- 0.3
eps1 <- 0.05
eps2 <- 0.05
# severity weight matrix
W <- matrix(c(0, 0, 1, 1.5, 2,
              0, 0, 5, 6, 7,
              0, 0, 2.5, 6, 8,
              0, 2, 3, 6, 9,
              0, 0, 1.5, 2, 2.5,
              0, 0, 0.5, 1, 1.5), nrow = 6, ncol = 5, byrow = TRUE)

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

# empty list to store results
# store tox response of different cohorts
tox_response <- list(d_1 = list(),
                     d_2 = list(),
                     d_3 = list(),
                     d_4 = list(),
                     d_5 = list(),
                     d_6 = list()) 
# store qTP of different cohorts
qTP <- list(d_1 = numeric(), 
            d_2 = numeric(), 
            d_3 = numeric(), 
            d_4 = numeric(),
            d_5 = numeric(),
            d_6 = numeric()) 
cohort_sequence <- numeric()
n_sequence <- numeric()
optimal_decision <- numeric()
optimal_decision_mtpi2 <- numeric()
# store DLT of different cohorts
DLT <- list(d_1 = numeric(),
            d_2 = numeric(),
            d_3 = numeric(),
            d_4 = numeric(),
            d_5 = numeric(),
            d_6 = numeric())
EqTP_record <- numeric()


# cohort 1 dose level 1 ----------------------------------------------------


## at dose level 1
# toxicity outcomes
tox_d1_1 <- matrix(c(0, 0, 1, 0, 0,
                     0, 0, 0, 0, 0,
                     0, 0, 1, 0, 0,
                     0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0,
                     0, 0, 1, 0, 0), nrow = 6, ncol = 5, byrow = TRUE)
tox_d1_2 <- matrix(c(0, 0, 1, 0, 0,
                     0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0), nrow = 6, ncol = 5, byrow = TRUE)
tox_d1_3 <- matrix(c(0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0,
                     0, 1, 0, 0, 0,
                     0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0), nrow = 6, ncol = 5, byrow = TRUE)
tox_d1 <- list(tox_d1_1, tox_d1_2, tox_d1_3)
tox_response[[1]] <- list.append(tox_response[[1]], tox_d1_1, tox_d1_2, tox_d1_3)
# calculate qTP
qTP_d1 <- numeric()
for(n in 1:length(tox_d1)){
  qTP_d1[n] <- sum(tox_d1[[n]] * W) / sum(apply(W, 1, max))
}
qTP[[1]] <- qTP_d1

# update dose-toxicity curve
qtp_sum<- unlist(sapply(qTP, sum))
qtp_n <- unlist(sapply(qTP, length))
l <- 1
data<-list(
  l = l,
  n = as.array(qtp_n[1:l]),
  sum_qtp = as.array(qtp_sum[1:l]),
  x = as.array(x[1:l]),
  a = -5,
  b = 5,
  c = 8
)
ab_MCMC <- stan('code/qtpi.stan', data = data, iter = 5000, refresh = 0)
params<-extract(ab_MCMC)
a_post <- params$alpha
b_post <- params$beta
# get decision
EqTP_i <- 1/(1+exp(-a_post-b_post*x[1]))
Mz_post <- as.numeric(table(cut(EqTP_i, int.bound)))
Mz_optimal <- which.max(Mz_post)
decision_optimal <- ifelse(Mz_optimal %in% qLI_index, "E", ifelse(Mz_optimal %in% qEI_index, "S", "D"))

# evaluate safety rule 1 and 2
EqTP_d1 <- 1/(1+exp(-a_post-b_post*x[1]))
mean(EqTP_d1 > qP_T) > 0.95 # safety 1 satisfied
EqTP_dnext <- 1/(1+exp(-a_post-b_post*x[2]))
mean(EqTP_dnext > qP_T) > 0.95 # safety 2 satisfied
# summary of dose level 1 (400 mg)
cohort_sequence[1] <- "400mg"
n_sequence[1] <- 3
optimal_decision[1] <- decision_optimal
EqTP_record[1] <- mean(EqTP_d1)

###############################################################################
# mTPI-2 for cohort 1, 400mg

d <- 1
mtpi2_result <- mtpi2(d, DLT, tox_response)
DLT <- mtpi2_result$DLT
Mz_post_mtpi2 <- mtpi2_result$Mz_post_mtpi2
Mz_optimal_mtpi2 <- mtpi2_result$Mz_optimal_mtpi2
decision_optimal_mtpi2 <- mtpi2_result$decision_optimal_mtpi2
optimal_decision_mtpi2[1] <- decision_optimal_mtpi2


# cohort 2 dose level 2 ----------------------------------------------------


## at dose level 2
# toxicity outcomes
tox_d2_1 <- matrix(c(0, 0, 1, 0, 0,
                     0, 0, 0, 0, 0,
                     0, 0, 1, 0, 0,
                     0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0,
                     0, 0, 1, 0, 0), nrow = 6, ncol = 5, byrow = TRUE)
tox_d2_2 <- matrix(c(0, 0, 0, 0, 0,
                     0, 0, 1, 0, 0,
                     0, 0, 1, 0, 0,
                     0, 1, 0, 0, 0,
                     0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0), nrow = 6, ncol = 5, byrow = TRUE)
tox_d2_3 <- matrix(c(0, 0, 1, 0, 0,
                     0, 0, 0, 0, 0,
                     0, 0, 1, 0, 0,
                     0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0,
                     0, 0, 1, 0, 0), nrow = 6, ncol = 5, byrow = TRUE)
tox_response[[2]] <- list.append(tox_response[[2]], tox_d2_1, tox_d2_2, tox_d2_3)
# calculate qTP
for(j in 1:(length(cohort_sequence)+1)){
  for(i in 1:length(tox_response[[j]])){
    qTP[[j]][i] <- sum(tox_response[[j]][[i]] * W) / sum(apply(W, 1, max))
  }
}
# update dose-toxicity curve
qtp_sum<- unlist(sapply(qTP, sum))
qtp_n <- unlist(sapply(qTP, length))
l <- 2
data<-list(
  l = l,
  n = as.array(qtp_n[1:l]),
  sum_qtp = as.array(qtp_sum[1:l]),
  x = as.array(x[1:l]),
  a = -5,
  b = 5,
  c = 8
)
ab_MCMC <- stan('code/qtpi.stan', data = data, iter = 5000, refresh = 0)
params<-extract(ab_MCMC)
a_post <- params$alpha
b_post <- params$beta
# get decision
EqTP_i <- 1/(1+exp(-a_post-b_post*x[2]))
Mz_post <- as.numeric(table(cut(EqTP_i, int.bound)))
Mz_optimal <- which.max(Mz_post)
decision_optimal <- ifelse(Mz_optimal %in% qLI_index, "E", ifelse(Mz_optimal %in% qEI_index, "S", "D"))
# evaluate safety rule 1 and 2
EqTP_d1 <- 1/(1+exp(-a_post-b_post*x[1]))
mean(EqTP_d1 > qP_T) > 0.95 # safety 1 satisfied
EqTP_dnext <- 1/(1+exp(-a_post-b_post*x[3]))
mean(EqTP_dnext > qP_T) > 0.95 # safety 2 satisfied
# summary of dose level 2 (500 mg)
cohort_sequence[2] <- "500mg"
n_sequence[2] <- 3
optimal_decision[2] <- decision_optimal
EqTP_d2 <- 1/(1+exp(-a_post-b_post*x[2]))
EqTP_record[2] <- mean(EqTP_d2)


###############################################################################
# mTPI-2 for cohort 2, 500mg

d <- 2
mtpi2_result <- mtpi2(d, DLT, tox_response)
DLT <- mtpi2_result$DLT
Mz_post_mtpi2 <- mtpi2_result$Mz_post_mtpi2
Mz_optimal_mtpi2 <- mtpi2_result$Mz_optimal_mtpi2
decision_optimal_mtpi2 <- mtpi2_result$decision_optimal_mtpi2
optimal_decision_mtpi2[2] <- decision_optimal_mtpi2



# cohort 3 dose level 3 ----------------------------------------------------


## at dose level 3
# toxicity outcomes
tox_d3_1 <- matrix(c(0, 0, 0, 0, 0,
                     0, 0, 1, 0, 0, 
                     0, 0, 1, 0, 0,
                     0, 0, 1, 0, 0, 
                     0, 0, 0, 0, 0, 
                     0, 0, 0, 0, 0), nrow = 6, ncol = 5, byrow = TRUE)
tox_d3_2 <- matrix(c(0, 0, 1, 0, 0,
                     0, 0, 0, 0, 0, 
                     0, 0, 0, 0, 0, 
                     0, 0, 1, 0, 0, 
                     0, 0, 1, 0, 0,
                     0, 0, 0, 0, 0), nrow = 6, ncol = 5, byrow = TRUE)
tox_d3_3 <- matrix(c(0, 0, 1, 0, 0,
                     0, 0, 0, 0, 0, 
                     0, 0, 1, 0, 0, 
                     0, 0, 0, 0, 0, 
                     0, 0, 0, 0, 0, 
                     0, 0, 1, 0, 0), nrow = 6, ncol = 5, byrow = TRUE)
tox_response[[3]] <- list.append(tox_response[[3]], tox_d3_1, tox_d3_2, tox_d3_3)
# calculate qTP
for(j in 1:(length(cohort_sequence)+1)){
  for(i in 1:length(tox_response[[j]])){
    qTP[[j]][i] <- sum(tox_response[[j]][[i]] * W) / sum(apply(W, 1, max))
  }
}
# update dose-toxicity curve
qtp_sum<- unlist(sapply(qTP, sum))
qtp_n <- unlist(sapply(qTP, length))
l <- 3
data<-list(
  l = l,
  n = as.array(qtp_n[1:l]),
  sum_qtp = as.array(qtp_sum[1:l]),
  x = as.array(x[1:l]),
  a = -5,
  b = 5,
  c = 8
)
ab_MCMC <- stan('code/qtpi.stan', data = data, iter = 5000, refresh = 0)
params<-extract(ab_MCMC)
a_post <- params$alpha
b_post <- params$beta
# get decision
EqTP_i <- 1/(1+exp(-a_post-b_post*x[3]))
Mz_post <- as.numeric(table(cut(EqTP_i, int.bound)))
Mz_optimal <- which.max(Mz_post)
decision_optimal <- ifelse(Mz_optimal %in% qLI_index, "E", ifelse(Mz_optimal %in% qEI_index, "S", "D"))
# evaluate safety rule 1 and 2
EqTP_d1 <- 1/(1+exp(-a_post-b_post*x[1]))
mean(EqTP_d1 > qP_T) > 0.95 # safety 1 satisfied
EqTP_dnext <- 1/(1+exp(-a_post-b_post*x[4]))
mean(EqTP_dnext > qP_T) > 0.95 # safety 2 satisfied
# summary of dose level 3 (600 mg)
cohort_sequence[3] <- "600mg"
n_sequence[3] <- 3
optimal_decision[3] <- decision_optimal
EqTP_d3 <- 1/(1+exp(-a_post-b_post*x[3]))
EqTP_record[3] <- mean(EqTP_d3)


################################################################################
# mTPI-2 for cohort 3, 600mg
d <- 3
mtpi2_result <- mtpi2(d, DLT, tox_response)
DLT <- mtpi2_result$DLT
Mz_post_mtpi2 <- mtpi2_result$Mz_post_mtpi2
Mz_optimal_mtpi2 <- mtpi2_result$Mz_optimal_mtpi2
decision_optimal_mtpi2 <- mtpi2_result$decision_optimal_mtpi2
optimal_decision_mtpi2[3] <- decision_optimal_mtpi2



# cohort 4 dose level 4 ----------------------------------------------------


## at dose level 4
# toxicity outcomes
tox_d4_1 <- matrix(c(0, 0, 0, 0, 0,
                     0, 0, 1, 0, 0, 
                     0, 0, 1, 0, 0, 
                     0, 0, 1, 0, 0, 
                     0, 0, 1, 0, 0, 
                     0, 0, 1, 0, 0), nrow = 6, ncol = 5, byrow = TRUE)
tox_d4_2 <- matrix(c(0, 0, 1, 0, 0,
                     0, 0, 0, 0, 0,
                     0, 0, 1, 0, 0,
                     0, 0, 1, 0, 0,
                     0, 0, 1, 0, 0,
                     0, 0, 0, 0, 0), nrow = 6, ncol = 5, byrow = TRUE)
tox_d4_3 <- matrix(c(0, 0, 1, 0, 0, 
                     0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0,
                     0, 0, 1, 0, 0, 
                     0, 0, 1, 0, 0,
                     0, 0, 1, 0, 0), nrow = 6, ncol = 5, byrow = TRUE)

tox_response[[4]] <- list.append(tox_response[[4]], tox_d4_1, tox_d4_2, tox_d4_3)
# calculate qTP
for(j in 1:(length(cohort_sequence)+1)){
  for(i in 1:length(tox_response[[j]])){
    qTP[[j]][i] <- sum(tox_response[[j]][[i]] * W) / sum(apply(W, 1, max))
  }
}
# update dose-toxicity curve
qtp_sum<- unlist(sapply(qTP, sum))
qtp_n <- unlist(sapply(qTP, length))
l <- 4
data<-list(
  l = l,
  n = as.array(qtp_n[1:l]),
  sum_qtp = as.array(qtp_sum[1:l]),
  x = as.array(x[1:l]),
  a = -5,
  b = 5,
  c = 8
)
ab_MCMC <- stan('code/qtpi.stan', data = data, iter = 5000, refresh = 0)
params<-extract(ab_MCMC)
a_post <- params$alpha
b_post <- params$beta
# get decision
EqTP_i <- 1/(1+exp(-a_post-b_post*x[4]))
Mz_post <- as.numeric(table(cut(EqTP_i, int.bound)))
Mz_optimal <- which.max(Mz_post)
decision_optimal <- ifelse(Mz_optimal %in% qLI_index, "E", ifelse(Mz_optimal %in% qEI_index, "S", "D"))
# evaluate safety rule 1 and 2
EqTP_d1 <- 1/(1+exp(-a_post-b_post*x[1]))
mean(EqTP_d1 > qP_T) > 0.95 # safety 1 satisfied
EqTP_dnext <- 1/(1+exp(-a_post-b_post*x[5]))
mean(EqTP_dnext > qP_T) > 0.95 # safety 2 satisfied
# summary of dose level 4 (700 mg)
cohort_sequence[4] <- "700mg"
n_sequence[4] <- 3
optimal_decision[4] <- decision_optimal
EqTP_d4 <- 1/(1+exp(-a_post-b_post*x[4]))
EqTP_record[4] <- mean(EqTP_d4)


################################################################################
# mTPI-2 for cohort 4, 700mg
d <- 4
mtpi2_result <- mtpi2(d, DLT, tox_response)
DLT <- mtpi2_result$DLT
Mz_post_mtpi2 <- mtpi2_result$Mz_post_mtpi2
Mz_optimal_mtpi2 <- mtpi2_result$Mz_optimal_mtpi2
decision_optimal_mtpi2 <- mtpi2_result$decision_optimal_mtpi2
optimal_decision_mtpi2[4] <- decision_optimal_mtpi2


# cohort 5 dose level 4 ----------------------------------------------------


## cohort 5, stay at dose level 4
tox_d4_1 <- matrix(c(0, 0, 0, 0, 0,
                     0, 0, 1, 0, 0, 
                     0, 0, 1, 0, 0, 
                     0, 0, 1, 0, 0, 
                     0, 0, 0, 0, 0, 
                     0, 1, 0, 0, 0), nrow = 6, ncol = 5, byrow = TRUE)
tox_d4_2 <- matrix(c(0, 0, 1, 0, 0,
                     0, 0, 0, 0, 0, 
                     0, 0, 0, 1, 0, 
                     0, 0, 0, 0, 0, 
                     0, 0, 1, 0, 0, 
                     0, 1, 0, 0, 0), nrow = 6, ncol = 5, byrow = TRUE)
tox_d4_3 <- matrix(c(0, 0, 1, 0, 0, 
                     0, 0, 0, 0, 0, 
                     0, 0, 1, 0, 0,
                     0, 0, 1, 0, 0, 
                     0, 0, 0, 0, 0,
                     0, 0, 1, 0, 0), nrow = 6, ncol = 5, byrow = TRUE)
tox_response[[4]] <- list.append(tox_response[[4]], tox_d4_1, tox_d4_2, tox_d4_3)
# calculate qTP
for(j in 1:4){
  for(i in 1:length(tox_response[[j]])){
    qTP[[j]][i] <- sum(tox_response[[j]][[i]] * W) / sum(apply(W, 1, max))
  }
}
# update dose-toxicity curve
qtp_sum<- unlist(sapply(qTP, sum))
qtp_n <- unlist(sapply(qTP, length))
l <- 4
data<-list(
  l = l,
  n = as.array(qtp_n[1:l]),
  sum_qtp = as.array(qtp_sum[1:l]),
  x = as.array(x[1:l]),
  a = -5,
  b = 5,
  c = 8
)
ab_MCMC <- stan('code/qtpi.stan', data = data, iter = 5000, refresh = 0)
params<-extract(ab_MCMC)
a_post <- params$alpha
b_post <- params$beta
# get decision
EqTP_i <- 1/(1+exp(-a_post-b_post*x[4]))
Mz_post <- as.numeric(table(cut(EqTP_i, int.bound)))
Mz_optimal <- which.max(Mz_post)
decision_optimal <- ifelse(Mz_optimal %in% qLI_index, "E", ifelse(Mz_optimal %in% qEI_index, "S", "D"))
# evaluate safety rule 1 and 2
EqTP_d1 <- 1/(1+exp(-a_post-b_post*x[1]))
mean(EqTP_d1 > qP_T) > 0.95 # safety 1 satisfied, no need to evaluate safety rule 2
# summary of dose level 4 (700 mg)
cohort_sequence[5] <- "700mg"
n_sequence[5] <- 3
optimal_decision[5] <- decision_optimal
EqTP_d4 <- 1/(1+exp(-a_post-b_post*x[4]))
EqTP_record[5] <- mean(EqTP_d4)

################################################################################
## at dose level 5
# toxicity outcomes
# tox_d5_1 <- matrix(c(0, 0, 0, 0, 0,
#                      0, 0, 1, 0, 0,
#                      0, 0, 1, 0, 0,
#                      0, 0, 1, 0, 0,
#                      0, 0, 0, 0, 0,
#                      0, 0, 1, 0, 0), nrow = 6, ncol = 5, byrow = TRUE)
# tox_d5_2 <- matrix(c(0, 0, 0, 1, 0,
#                      0, 0, 0, 0, 0,
#                      0, 0, 0, 0, 1,
#                      0, 0, 0, 1, 0,
#                      0, 0, 0, 1, 0,
#                      0, 0, 1, 0, 0), nrow = 6, ncol = 5, byrow = TRUE)
# tox_d5_3 <- matrix(c(0, 0, 1, 0, 0,
#                      0, 0, 0, 0, 0,
#                      0, 0, 1, 0, 0,
#                      0, 0, 1, 0, 0, 
#                      0, 0, 0, 0, 0, 
#                      0, 0, 1, 0, 0), nrow = 6, ncol = 5, byrow = TRUE)
# tox_response[[5]] <- list.append(tox_response[[5]], tox_d5_1, tox_d5_2, tox_d5_3)
# # calculate qTP
# for(j in 1:(length(cohort_sequence)+1)){
#   for(i in 1:length(tox_response[[j]])){
#     qTP[[j]][i] <- sum(tox_response[[j]][[i]] * W) / sum(apply(W, 1, max))
#   }
# }
# # update dose-toxicity curve
# ab_MCMC <- qTPI_MCMC(qTP = qTP, x = x, iter=10000)
# a_post <- ab_MCMC$a_posterior
# b_post <- ab_MCMC$b_posterior
# # get decision
# Mz_post <- rep(0, length(int.bound)-1)
# for (i in 1:length(b_post)){
#   a_i <- a_post[i]
#   b_i <- b_post[i]
#   EqTP_i <- 1/(1+exp(-a_i-b_i*x[5]))
#   for (index in 1:(length(int.bound) - 1)) {
#     low_bound <- int.bound[index]
#     up_bound <- int.bound[index + 1]
#     if (EqTP_i >= low_bound & EqTP_i < up_bound){
#       Mz_post[index] <- Mz_post[index] + 1
#     }
#   }
# }
# Mz_optimal <- which.max(Mz_post)
# decision_optimal <- ifelse(Mz_optimal %in% qLI_index, "E", ifelse(Mz_optimal %in% qEI_index, "S", "D"))
# # evaluate safety rule 1 and 2
# EqTP_d1 <- 1/(1+exp(-a_post-b_post*x[1]))
# mean(EqTP_d1 > qP_T) > 0.95 # safety 1 satisfied, no need to evaluate safety rule 2
# # summary of dose level 4 (700 mg)
# cohort_sequence[5] <- "800mg"
# n_sequence[5] <- 3
# optimal_decision[5] <- decision_optimal
# EqTP_d5 <- 1/(1+exp(-a_post-b_post*x[5]))
# EqTP_record[5] <- mean(EqTP_d5)


################################################################################
# mTPI-2 for cohort 5, 800mg
# d <- 5
# mtpi2_result <- mtpi2(d, DLT, tox_response)
# DLT <- mtpi2_result$DLT
# Mz_post_mtpi2 <- mtpi2_result$Mz_post_mtpi2
# Mz_optimal_mtpi2 <- mtpi2_result$Mz_optimal_mtpi2
# decision_optimal_mtpi2 <- mtpi2_result$decision_optimal_mtpi2
# optimal_decision_mtpi2[5] <- decision_optimal_mtpi2


# cohort 6 dose level 4 ---------------------------------------------------


## cohort 6, stay at dose level 4
# toxicity outcomes
tox_d4_1 <- matrix(c(0, 0, 1, 0, 0,
                     0, 0, 0, 0, 0, 
                     0, 0, 1, 0, 0, 
                     0, 0, 1, 0, 0, 
                     0, 0, 1, 0, 0, 
                     0, 0, 0, 1, 0), nrow = 6, ncol = 5, byrow = TRUE)
tox_d4_2 <- matrix(c(0, 0, 1, 0, 0,
                     0, 0, 0, 0, 0, 
                     0, 0, 1, 0, 0, 
                     0, 0, 1, 0, 0, 
                     0, 0, 1, 0, 0, 
                     0, 0, 0, 0, 0), nrow = 6, ncol = 5, byrow = TRUE)
tox_d4_3 <- matrix(c(0, 0, 0, 1, 0, 
                     0, 0, 0, 0, 0, 
                     0, 0, 1, 0, 0,
                     0, 0, 1, 0, 0, 
                     0, 0, 0, 0, 0,
                     0, 0, 1, 0, 0), nrow = 6, ncol = 5, byrow = TRUE)
tox_response[[4]] <- list.append(tox_response[[4]], tox_d4_1, tox_d4_2, tox_d4_3)
# calculate qTP
for(j in 1:4){
  for(i in 1:length(tox_response[[j]])){
    qTP[[j]][i] <- sum(tox_response[[j]][[i]] * W) / sum(apply(W, 1, max))
  }
}
# update dose-toxicity curve
qtp_sum<- unlist(sapply(qTP, sum))
qtp_n <- unlist(sapply(qTP, length))
l <- 4
data<-list(
  l = l,
  n = as.array(qtp_n[1:l]),
  sum_qtp = as.array(qtp_sum[1:l]),
  x = as.array(x[1:l]),
  a = -5,
  b = 5,
  c = 8
)
ab_MCMC <- stan('code/qtpi.stan', data = data, iter = 5000, refresh = 0)
params<-extract(ab_MCMC)
a_post <- params$alpha
b_post <- params$beta
# get decision
EqTP_i <- 1/(1+exp(-a_post-b_post*x[4]))
Mz_post <- as.numeric(table(cut(EqTP_i, int.bound)))
Mz_optimal <- which.max(Mz_post)
decision_optimal <- ifelse(Mz_optimal %in% qLI_index, "E", ifelse(Mz_optimal %in% qEI_index, "S", "D"))
# evaluate safety rule 1 and 2
EqTP_d1 <- 1/(1+exp(-a_post-b_post*x[1]))
mean(EqTP_d1 > qP_T) > 0.95 # safety 1 satisfied, no need to evaluate safety rule 2
# summary of dose level 4 (700 mg)
cohort_sequence[6] <- "700mg"
n_sequence[6] <- 3
optimal_decision[6] <- decision_optimal
EqTP_d4 <- 1/(1+exp(-a_post-b_post*x[4]))
EqTP_record[6] <- mean(EqTP_d4)


################################################################################
# mTPI-2 for cohort 6, 700mg
d <- 4
mtpi2_result <- mtpi2(d, DLT, tox_response)
DLT <- mtpi2_result$DLT
Mz_post_mtpi2 <- mtpi2_result$Mz_post_mtpi2
Mz_optimal_mtpi2 <- mtpi2_result$Mz_optimal_mtpi2
decision_optimal_mtpi2 <- mtpi2_result$decision_optimal_mtpi2
optimal_decision_mtpi2[6] <- decision_optimal_mtpi2


# cohort 7 dose level 4 ---------------------------------------------------


## cohort 7,  Stay at dose level 4
# toxicity outcomes
tox_d4_1 <- matrix(c(0, 0, 1, 0, 0, 
                     0, 0, 0, 0, 0, 
                     0, 0, 0, 0, 0, 
                     0, 1, 0, 0, 0, 
                     0, 0, 1, 0, 0, 
                     0, 0, 1, 0, 0), nrow = 6, ncol = 5, byrow = TRUE)
tox_d4_2 <- matrix(c(0, 0, 0, 0, 0,
                     0, 0, 1, 0, 0, 
                     0, 0, 1, 0, 0, 
                     0, 0, 0, 0, 0, 
                     0, 0, 0, 0, 0,
                     0, 0, 1, 0, 0), nrow = 6, ncol = 5, byrow = TRUE)
tox_d4_3 <- matrix(c(0, 0, 0, 1, 0,
                     0, 0, 0, 0, 0, 
                     0, 0, 1, 0, 0, 
                     0, 0, 0, 1, 0, 
                     0, 0, 1, 0, 0,
                     0, 0, 0, 1, 0), nrow = 6, ncol = 5, byrow = TRUE)
tox_response[[4]] <- list.append(tox_response[[4]], tox_d4_1, tox_d4_2, tox_d4_3)
# calculate qTP
for(j in 1:4){
  for(i in 1:length(tox_response[[j]])){
    qTP[[j]][i] <- sum(tox_response[[j]][[i]] * W) / sum(apply(W, 1, max))
  }
}
# update dose-toxicity curve
qtp_sum<- unlist(sapply(qTP, sum))
qtp_n <- unlist(sapply(qTP, length))
l <- 4
data<-list(
  l = l,
  n = as.array(qtp_n[1:l]),
  sum_qtp = as.array(qtp_sum[1:l]),
  x = as.array(x[1:l]),
  a = -5,
  b = 5,
  c = 8
)
ab_MCMC <- stan('code/qtpi.stan', data = data, iter = 5000, refresh = 0)
params<-extract(ab_MCMC)
a_post <- params$alpha
b_post <- params$beta
# get decision
EqTP_i <- 1/(1+exp(-a_post-b_post*x[4]))
Mz_post <- as.numeric(table(cut(EqTP_i, int.bound)))
Mz_optimal <- which.max(Mz_post)
decision_optimal <- ifelse(Mz_optimal %in% qLI_index, "E", ifelse(Mz_optimal %in% qEI_index, "S", "D"))
# evaluate safety rule 1 and 2
EqTP_d1 <- 1/(1+exp(-a_post-b_post*x[1]))
mean(EqTP_d1 > qP_T) > 0.95 # safety 1 satisfied, no need to evaluate safety rule 2
# summary of dose level 4 (700 mg)
cohort_sequence[7] <- "700mg"
n_sequence[7] <- 3
optimal_decision[7] <- decision_optimal
EqTP_d4 <- 1/(1+exp(-a_post-b_post*x[4]))
EqTP_record[7] <- mean(EqTP_d4)


################################################################################
# mTPI-2 for cohort 7, 700mg
d <- 4
mtpi2_result <- mtpi2(d, DLT, tox_response)
DLT <- mtpi2_result$DLT
Mz_post_mtpi2 <- mtpi2_result$Mz_post_mtpi2
Mz_optimal_mtpi2 <- mtpi2_result$Mz_optimal_mtpi2
decision_optimal_mtpi2 <- mtpi2_result$decision_optimal_mtpi2
optimal_decision_mtpi2[7] <- decision_optimal_mtpi2


# cohort 8 dose level 4 ---------------------------------------------------


## cohort 8, Stay at dose level 4 again
# toxicity outcomes
tox_d4_1 <- matrix(c(0, 0, 0, 0, 0,
                     0, 0, 1, 0, 0, 
                     0, 0, 0, 0, 0, 
                     0, 1, 0, 0, 0,
                     0, 0, 0, 0, 0, 
                     0, 0, 1, 0, 0), nrow = 6, ncol = 5, byrow = TRUE)
tox_d4_2 <- matrix(c(0, 0, 0, 1, 0,
                     0, 0, 0, 0, 0,
                     0, 0, 1, 0, 0,
                     0, 0, 0, 0, 0,
                     0, 0, 1, 0, 0,
                     0, 0, 1, 0, 0), nrow = 6, ncol = 5, byrow = TRUE)
tox_d4_3 <- matrix(c(0, 0, 1, 0, 0,
                     0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0,
                     0, 0, 1, 0, 0, 
                     0, 0, 1, 0, 0,
                     0, 0, 1, 0, 0), nrow = 6, ncol = 5, byrow = TRUE)
tox_response[[4]] <- list.append(tox_response[[4]], tox_d4_1, tox_d4_2, tox_d4_3)
# calculate qTP
for(j in 1:4){
  for(i in 1:length(tox_response[[j]])){
    qTP[[j]][i] <- sum(tox_response[[j]][[i]] * W) / sum(apply(W, 1, max))
  }
}
# update dose-toxicity curve
qtp_sum<- unlist(sapply(qTP, sum))
qtp_n <- unlist(sapply(qTP, length))
l <- 4
data<-list(
  l = l,
  n = as.array(qtp_n[1:l]),
  sum_qtp = as.array(qtp_sum[1:l]),
  x = as.array(x[1:l]),
  a = -5,
  b = 5,
  c = 8
)
ab_MCMC <- stan('code/qtpi.stan', data = data, iter = 5000, refresh = 0)
params<-extract(ab_MCMC)
a_post <- params$alpha
b_post <- params$beta
# get decision
EqTP_i <- 1/(1+exp(-a_post-b_post*x[4]))
Mz_post <- as.numeric(table(cut(EqTP_i, int.bound)))
Mz_optimal <- which.max(Mz_post)
decision_optimal <- ifelse(Mz_optimal %in% qLI_index, "E", ifelse(Mz_optimal %in% qEI_index, "S", "D"))
# evaluate safety rule 1 and 2
EqTP_d1 <- 1/(1+exp(-a_post-b_post*x[1]))
mean(EqTP_d1 > qP_T) > 0.95 # safety 1 satisfied, no need to evaluate safety rule 2
# summary of dose level 4 (700 mg)
cohort_sequence[8] <- "700mg"
n_sequence[8] <- 3
optimal_decision[8] <- decision_optimal
EqTP_d4 <- 1/(1+exp(-a_post-b_post*x[4]))
EqTP_record[8] <- mean(EqTP_d4)


################################################################################
# mTPI-2 for cohort 8, 700mg
d <- 4
mtpi2_result <- mtpi2(d, DLT, tox_response)
DLT <- mtpi2_result$DLT
Mz_post_mtpi2 <- mtpi2_result$Mz_post_mtpi2
Mz_optimal_mtpi2 <- mtpi2_result$Mz_optimal_mtpi2
decision_optimal_mtpi2 <- mtpi2_result$decision_optimal_mtpi2
optimal_decision_mtpi2[8] <- decision_optimal_mtpi2


# cohort 9 dose level 4 ---------------------------------------------------


## cohort 9, Stay at dose level 4 again
# toxicity outcomes
tox_d4_1 <- matrix(c(0, 0, 0, 0, 0,
                     0, 0, 1, 0, 0,
                     0, 0, 0, 0, 0,
                     0, 0, 1, 0, 0,
                     0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0), nrow = 6, ncol = 5, byrow = TRUE)
tox_d4_2 <- matrix(c(0, 0, 1, 0, 0, 
                     0, 0, 0, 0, 0, 
                     0, 0, 1, 0, 0, 
                     0, 1, 0, 0, 0, 
                     0, 0, 0, 0, 0, 
                     0, 0, 1, 0, 0), nrow = 6, ncol = 5, byrow = TRUE)
tox_d4_3 <- matrix(c(0, 0, 1, 0, 0,
                     0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 
                     0, 1, 0, 0, 0,
                     0, 0, 0, 1, 0, 
                     0, 0, 0, 1, 0), nrow = 6, ncol = 5, byrow = TRUE)
tox_response[[4]] <- list.append(tox_response[[4]], tox_d4_1, tox_d4_2, tox_d4_3)
# calculate qTP
for(j in 1:4){
  for(i in 1:length(tox_response[[j]])){
    qTP[[j]][i] <- sum(tox_response[[j]][[i]] * W) / sum(apply(W, 1, max))
  }
}
# update dose-toxicity curve
qtp_sum<- unlist(sapply(qTP, sum))
qtp_n <- unlist(sapply(qTP, length))
l <- 4
data<-list(
  l = l,
  n = as.array(qtp_n[1:l]),
  sum_qtp = as.array(qtp_sum[1:l]),
  x = as.array(x[1:l]),
  a = -5,
  b = 5,
  c = 8
)
ab_MCMC <- stan('code/qtpi.stan', data = data, iter = 5000, refresh = 0)
params<-extract(ab_MCMC)
a_post <- params$alpha
b_post <- params$beta
# get decision
EqTP_i <- 1/(1+exp(-a_post-b_post*x[4]))
Mz_post <- as.numeric(table(cut(EqTP_i, int.bound)))
Mz_optimal <- which.max(Mz_post)
decision_optimal <- ifelse(Mz_optimal %in% qLI_index, "E", ifelse(Mz_optimal %in% qEI_index, "S", "D"))
# evaluate safety rule 1 and 2
EqTP_d1 <- 1/(1+exp(-a_post-b_post*x[1]))
mean(EqTP_d1 > qP_T) > 0.95 # safety 1 satisfied, no need to evaluate safety rule 2
# summary of dose level 4 (700 mg)
cohort_sequence[9] <- "700mg"
n_sequence[9] <- 3
optimal_decision[9] <- decision_optimal
EqTP_d4 <- 1/(1+exp(-a_post-b_post*x[4]))
EqTP_record[9] <- mean(EqTP_d4)


################################################################################
# mTPI-2 for cohort 9, 700mg
d <- 4
mtpi2_result <- mtpi2(d, DLT, tox_response)
DLT <- mtpi2_result$DLT
Mz_post_mtpi2 <- mtpi2_result$Mz_post_mtpi2
Mz_optimal_mtpi2 <- mtpi2_result$Mz_optimal_mtpi2
decision_optimal_mtpi2 <- mtpi2_result$decision_optimal_mtpi2
optimal_decision_mtpi2[9] <- decision_optimal_mtpi2

################################################################################
## After trial finish
# select MTD (700)
# percent treated at MTD :15/27
# percent treated above MTD: 3/27



