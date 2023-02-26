################################################
# Dosing Scenarios for Simulation
################################################
# function for generating probability of grade0to4 across all 6 levels (for one tox type?)
p_all_doses <- function(p_initial, D, trans_matrix){
  p_all <- matrix(0, nrow=D, ncol=length(p_initial))
  p_all[1, ] <- p_initial
  for(r in 2:D){
    for(c in 1:length(p_initial)){
      p_all[r,c] <- p_all[r-1, ] %*% trans_matrix[,c]
    }
  }
  return(p_all)
}

# function to find the MTD based on qTP
MTD_qTP <- function(p_renal, p_neurological, p_hematological, D, W, P_T, eps2){
  W_max <- sum(apply(W, 1, max)) # sum of the maximum of each row of W
  qTP <- numeric()
  for(d in 1:D){
    dose_d <- rbind(p_renal[d,],
                    p_neurological[d,],
                    p_hematological[d,])
    qTP[d] <- sum(dose_d * W) / W_max
  }
  if(sum(qTP > P_T+eps2) == D){ # qTP's of all doses greater than 0.35
    MTD <- 0 # no MTD
  }else{
    qTP_candidate <- qTP[qTP <= P_T+eps2] # for doses whose qTP <= 0.35
    MTD <- which.min(abs(qTP_candidate - P_T)) # choose the dose whose qTP is closest to 0.3
  }
  return(list(qTP = qTP, MTD_qTP = MTD))
}

# function to find MTD based on binary p
MTD_p <- function(p_renal, p_neurological, p_hematological, D, P_T, eps2){
  binary_renal <- numeric()
  binary_neurological <- numeric()
  binary_hematological <- numeric()
  binary_p <- numeric()
  for(d in 1:D){
    binary_renal[d] <- p_renal[d,4] + p_renal[d,5]
    binary_neurological[d] <- p_neurological[d,4] + p_neurological[d,5] # grade >=3 as DLT
    binary_hematological[d] <- p_hematological[d,5] # grade >=4 as DLT
    binary_p[d] <- 1-(1-binary_renal[d])*(1-binary_neurological[d])*(1-binary_hematological[d])
  }
  # determine MTD by binary p
  if(sum(binary_p > P_T+eps2) == D){
    MTD_p <- 0
  }else{
    binary_p_candidate <- binary_p[binary_p <= P_T+eps2]
    MTD_p <- which.min(abs(binary_p_candidate - P_T))
  }
  return(list(p=binary_p,
              MTD=MTD_p))
}

################################################
# Generating Dosing Scenarios for Simulation
################################################
# there are a total of 27 dosing scenarios 
# for 3 toxicity types with (L, M, H) levels of transition matrix
Gen_Scn <- function(transition_level = c("low", "low", "low"),
                    p_r, p_n, p_h, D, W, P_T, eps2){
  # for renal tox
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
  if(transition_level[1] == "low"){
    p_renal <- p_all_doses(p_initial = p_r, D = D, trans_matrix = trans_matrix_L)
  }else if(transition_level[1] == "moderate"){
    p_renal <- p_all_doses(p_initial = p_r, D = D, trans_matrix = trans_matrix_M)
  }else if(transition_level[1] == "high"){
    p_renal <- p_all_doses(p_initial = p_r, D = D, trans_matrix = trans_matrix_H)
  }
  # for neurological
  if(transition_level[2] == "low"){
    p_neurological <- p_all_doses(p_initial = p_n, D = D, trans_matrix = trans_matrix_L)
  }else if(transition_level[2] == "moderate"){
    p_neurological <- p_all_doses(p_initial = p_n, D = D, trans_matrix = trans_matrix_M)
  }else if(transition_level[2] == "high"){
    p_neurological <- p_all_doses(p_initial = p_n, D = D, trans_matrix = trans_matrix_H)
  }
  # for hematological
  if(transition_level[3] == "low"){
    p_hematological <- p_all_doses(p_initial = p_h, D = D, trans_matrix = trans_matrix_L)
  }else if(transition_level[3] == "moderate"){
    p_hematological <- p_all_doses(p_initial = p_h, D = D, trans_matrix = trans_matrix_M)
  }else if(transition_level[3] == "high"){
    p_hematological <- p_all_doses(p_initial = p_h, D = D, trans_matrix = trans_matrix_H)
  }
  # calculate MTD based on qTP
  MTD_by_qTP <- MTD_qTP(p_renal,p_neurological,p_hematological,D,W,P_T,eps2)
  MTD_by_p <- MTD_p(p_renal, p_neurological, p_hematological, D, P_T, eps2)
  return(list(MTD_by_qTP = MTD_by_qTP$MTD,
              qTP = MTD_by_qTP$qTP,
              MTD_by_p = MTD_by_p$MTD,
              p = MTD_by_p$p,
              p_all = rbind(p_renal,p_neurological,p_hematological)))
}

###################################################################
# record qTP, p, MTD(selected by qTP), MTD(selected by p)
###################################################################
record_qTP <- function(p_r, p_n, p_h, W, P_T, eps2){
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
  
  qTP <- matrix(numeric(), nrow = 27, ncol = 6)
  p <- matrix(numeric(), nrow = 27, ncol = 6)
  MTD_by_qTP <- numeric(27)
  MTD_by_p <- numeric(27)
  
  for ( i in 1: length(scenario_list)){
    transition_level <- scenario_list[[i]]
    Scn <- Gen_Scn(transition_level = transition_level, p_r = p_r, p_n = p_n, p_h = p_h, D = 6, W = W, P_T, eps2)
    qTP[i,] <- Scn$qTP
    p[i,] <- Scn$p
    MTD_by_qTP[i] <- Scn$MTD_by_qTP
    MTD_by_p[i] <- Scn$MTD_by_p
  }

  return(list(
    qTP = qTP,
    p = p,
    MTD_by_qTP = MTD_by_qTP,
    MTD_by_p = MTD_by_p
  ))
}

