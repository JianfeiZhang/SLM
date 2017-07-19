###################################################################
# SLM.R
# Copyright (C) 2017 Université de Sherbrooke. All Rights Reserved
# This file is a part of CHUS clinical project -
# survival learning machine
# Written by Jianfei Zhang <jianfei.zhang@live.ca>
# May 2017
###################################################################

library(caret)
library(survival)
library(GA)
library(pec)

##################
# DATA READ
##################

data <-
  read.table(Data.path,
             header = TRUE,
             sep = ";")
# 38 risk factors
feat <- c(9, 15:51)
event <- data[, 4]

# stratified 5-fold cross-validation
folds = 5
cvfold <-
  createFolds(factor(event), folds, list = TRUE, returnTrain = TRUE)

AUC <- 0
CI <- 0
MSE <- 0

start.time <- Sys.time()

for (k in 1:folds) {
  train_data <- data[cvfold[[k]], feat]
  test_data <- data[-cvfold[[k]], feat]
  # follow-up time
  train_time <- data[cvfold[[k]], 5]
  test_time <- data[-cvfold[[k]], 5]
  # survival status
  train_event <- data[cvfold[[k]], 4]
  test_event <- data[-cvfold[[k]], 4]
  
  ##############################
  # COMPUTE THE COX MODEL
  ##############################
  surv <- Surv(train_time, train_event == 1)
  model.cox <-
    coxph(surv ~ ., data = data.frame(train_data), ties = "exact")
  beta <- as.vector(model.cox$coefficients)
  
  
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
    f - penalty1 - penalty2   # fitness function value
  }
  
  res <- ga(
    "real-valued",
    fitness = fitness,
    min = lb,
    max = ub,
    maxiter = 1000
  )
  
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

end.time <- Sys.time()

print(round(
  c(AUC / folds, CI / folds, MSE / folds, end.time - start.time),
  digits = 3
))
