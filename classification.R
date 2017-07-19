###################################################################
# classification.R
# Copyright (C) 2017 Université de Sherbrooke. All Rights Reserved
# This file is a part of CHUS clinical project -
# survival learning machine
# Written by Jianfei Zhang <jianfei.zhang@live.ca>
# May 2017
###################################################################

library(caret)
library(e1071)

##################
# DATA READ
##################

data <-
  read.table(Data.path,
             header = TRUE,
             sep = ";")
# the full set or subset of selected risk factors
feat <- c(9, 15:51)

event <- data[, 4]

# stratified 5-fold cross-validation
folds = 5
cvfold <-
  createFolds(factor(event), folds, list = TRUE, returnTrain = TRUE)

recall <- 0
spe <- 0
pre <- 0

for (k in 1:folds) {
  train_data <- data[cvfold[[k]], feat]
  test_data <- data[-cvfold[[k]], feat]
  # follow-up time
  train_time <- data[cvfold[[k]], 5]
  test_time <- data[-cvfold[[k]], 5]
  # survival status
  train_event <- data[cvfold[[k]], 4]
  test_event <- data[-cvfold[[k]], 4]
  
  
  ###############################################
  # Classification Model
  ###############################################
  # glm model
  model <-
    glm(as.factor(train_event) ~ . ,
        family = binomial(logit),
        data = train_data)
  # svm model
  # model<-
  #   svm(as.factor(train_event) ~ ., data = train_data, kernel = "linear")
  
  # glm classification
  pred <-
    ifelse(predict(model, newdata = test_data) > 0, '1', '0')
  
  # svm classification
  # pred <- predict(model, newdata = test_data)
  table(pred, test_event)
  
  # sensitivity (recall)
  recall <-
    recall + sensitivity(as.factor(pred), as.factor(test_event), positive = 1)
  
  # specificity
  spe <-
    spe + specificity(as.factor(pred), as.factor(test_event), negative = 0)
  
  # precision
  pre <-
    pre + posPredValue(as.factor(pred), as.factor(test_event), positive = 1)
  
}

# average 5CV recall, specificity, precision
getAvg <- function(metric, num_cv) {
  return(round(metric / num_cv))
}

recall = recall / folds
spe = spe / folds
pre = pre / folds

# average 5CV F1-score
F1 <- (2 * pre * recall) / (pre + recall)

print(round(c(recall, spe, F1), digits = 3))
