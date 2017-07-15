library(caret)
library(survival)
library(GA)
library(doParallel)
library(foreach)
library(iterators)
library(parallel)
library(pec)

##################
# DATA READ
##################
data <-
  read.table(
    "/home/jianfeizhang/Documents/rstudio/SelectFactor/data.csv",
    header = TRUE,
    sep = ";"
  )
# 38 risk factors
feat <- c(9, 15:51)
event <- data[, 4]

# 5-fold cross-validation
folds = 5
cvfold <-
  createFolds(factor(event), folds, list = TRUE, returnTrain = TRUE)

AUC <- 0
CI <- 0
MSE <- 0

for (k in 1:folds) {
  train_data <- data[cvfold[[k]], feat]
  test_data <- data[-cvfold[[k]], feat]
  # follow-up time
  train_time <- data[cvfold[[k]], 5]
  test_time <- data[-cvfold[[k]], 5]
  # survival status
  train_event <- data[cvfold[[k]], 4]
  test_event <- data[-cvfold[[k]], 4]
  #################################################
  # Build the failure and failure-free sets
  #################################################
  # divide_set <- function(train_event) {
  #   E0 <- c()
  #   E1 <- c()
  #   for (i in 1:length(train_event)) {
  #     if (train_event[i] == 1) {
  #       # find indices
  #       E1 <- c(E1 , cvfold[[1]][i])
  #     } else{
  #       E0 <- c(E0 , cvfold[[1]][i])
  #     }
  #   }
  #   ans <- list("E0" = E0,
  #               "E1" = E1)
  #   return(ans)
  # }
  #
  # tools <- divide_set(train_event)
  
  ##############################
  # COMPUTE THE COX MODEL
  ##############################
  surv <- Surv(train_time, train_event == 1)
  model.cox <-
    coxph(surv ~ ., data = data.frame(train_data), ties = "exact")
  beta <- as.vector(model.cox$coefficients)
  
  
  ####################################################
  # FUNCTIONS MODELIZING OUR OPTIMIZATION PROBLEM
  ####################################################
  
  # weighted value of link function
  # survival_p <- function(w , train_data) {
  #   E0 <- c()
  #   E1 <- c()
  #   for (i in cvfold[[1]]) {
  #     # event=1
  #     if (data[i,4] == 1) {
  #       E1 <- c(E1 , sum(beta * getValues(i, train_data) * w))
  #     } else{
  #       E0 <- c(E0 , sum(beta * getValues(i, train_data) * w))
  #     }
  #   }
  #   ans <- list("E0" = E0,
  #               "E1" = E1)
  #   return(ans)
  # }
  
  # sigmoid constant
  theta <- 1e6
  
  # Dirichlet hyperparameter
  alpha <- 1.5
  
  objective_function <- function(w) {
    ans <- 0
    for (i in 1:length(train_data[, 1])) {
      # failure patients
      if (train_event[i] == 1) {
        bdw_i <- sum(beta * as.numeric(train_data[i, ]) * w)
        for (j in 1:length(train_data[, 1])) {
          if (train_time[j] > train_time[i]) {
            bdw_j <- sum(beta * as.numeric(train_data[j, ]) * w)
            sigm <-
              1 / (1 + exp(theta * ((bdw_j - bdw_i) - 1 / sqrt(theta)
              )))
            ans <- ans + sum(sigm)
          }
        }
      }
    }
    b <- length(feat) * gamma(alpha) / gamma(length(feat) * alpha)
    log_Dir <- -log(b) + (alpha - 1) * sum(log(w))
    return(ans + log_Dir)
  }
  
  # grad_objective_function <- function( w ){
  #   ans <- rep(0,length(w))
  #   npair <- ( length(tools$E0) * length(tools$E1) )
  #   for( i in tools$E0 ){
  #     xi <- getPatient(i,datax)
  #     for( j in tools$E1 ){
  #       xj <- getPatient(j,datax)
  #       delta_ij <- xi - xj
  #       bdw <- sum( beta * delta_ij * w )
  #       exptheta <- exp( theta * ( bdw - 1 / sqrt(theta) ) )
  #       #print(exptheta)
  #       ans <- ans + theta * exptheta * ( beta * delta_ij ) / ( 1 + exptheta )^2
  #       #print(ans)
  #     }
  #   }
  #   ans <- ans / npair
  #   return(ans)
  # }
  
  ################################################
  # FUNCTIONS TO SOLVE THE OPTIMIZATION PROBLEM
  ################################################
  
  # the function we want to minimize
  eval_f <- function(w) {
    objective_function(w)
  }
  
  # tolrance
  mu <- 0.02
  
  # the equality constraint
  eval_g_eq1 <- function(w) {
    sum(w) - length(w) - mu
  }
  
  eval_g_eq2 <- function(w) {
    length(w) - sum(w) - mu
  }
  
  # initial values
  # w <- rep(0.5, length(feat))
  
  # lower and upper bounds
  lb <- as.numeric(rep(0, length(feat)))
  ub <- as.numeric(rep(1, length(feat)))
  
  fitness <- function(w)
  {
    f <- -objective_function(w)  # maximise -f(x)
    print(c(f, sum(w)))
    # penalty term
    penalty1 <-
      max(eval_g_eq1(w), 0)   # penalisation for 1st inequality constraint
    penalty2 <-
      max(eval_g_eq2(w), 0)   # penalisation for 2nd inequality constraint
    f - penalty1 - penalty2            # fitness function value
  }
  
  res <- ga(
    "real-valued",
    fitness = fitness,
    min = lb,
    max = ub,
    maxiter = 1000
  )
  
  # library(nloptr)
  # res<-nloptr(w0,eval_f = eval_f,lb = lb,ub = ub)
  
  # objective_function(rep(1,39))*57400
  # objective_function(w)*57400
  # write.table(w,file = "w.txt",sep = ";")
  
  ###############################################
  # Model Evaluation in terms of AUC and C-index
  ###############################################
  getAUCCI <- function(w) {
    pair_AUC <- 0
    correct_pair_AUC <- 0
    pair_CI <- 0
    correct_pair_CI <- 0
    for (i in 1:length(test_data[, 1])) {
      # for given failure patients
      if (test_event[i] == 1) {
        bdw_i <- sum(beta * as.numeric(test_data[i, ]) * w)
        for (j in 1:length(test_data[, 1])) {
          bdw_j <- sum(beta * as.numeric(test_data[j, ]) * w)
          
          #compute AUC
          if (test_event[j] == 0) {
            pair_AUC = pair_AUC + 1
            if (bdw_j > bdw_i) {
              correct_pair_AUC = correct_pair_AUC + 1
            }
          }
          
          # compute C-index
          if (test_time[j] > test_time[i]) {
            pair_CI = pair_CI + 1
            if (bdw_j > bdw_i) {
              correct_pair_CI = correct_pair_CI + 1
            }
          }
        }
      }
    }
    return(c(correct_pair_AUC / pair_AUC, correct_pair_CI / pair_CI))
  }
  
  ###############################################
  # Model Evaluation in terms of MSE
  ###############################################
  getMSE <- function(w) {
    loss <- 0
    for (i in 1:length(test_data[, 1])) {
      presurv <-
        predictSurvProb(model.cox, newdata = test_data[i,], times = seq(365))
      loss = loss + (1 - as.numeric(test_event[i]) - presurv) ^ 2
    }
    return(loss)
  }
  
  AUC <- getAUCCI(res@solution)[1]
  CI <- getAUCCI(res@solution)[2]
  MSE <- MSE + getMSE(res@solution)
  
}
print(c(AUC / 5, CI / 5, MSE / 5))
