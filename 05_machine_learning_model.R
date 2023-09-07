###h2o package for feature selection
library(h2o)
features <- colnames(s.data)
s.data<-s.data[,marker]
h2o.init()
train.data.h2o<-as.h2o(s.data)

anomaly_model <- h2o.deeplearning(x = features, 
                                  training_frame = train.data.h2o,
                                  activation = "Tanh", 
                                  autoencoder = TRUE, 
                                  hidden = 100,
                                  epochs = 100,
                                  variable_importances=TRUE,
                                  seed=20230731,
                                  standardize=TRUE)

head(anomaly_model@model$variable_importances)

import_mat2<-anomaly_model@model$variable_importances

###MLP model by neuralnet
mlpcla <- neuralnet(class~., 
                    data = x_train,
                    stepmax = 500000,
                    rep = 3,
                    hidden = c(8), ## 隐藏层神经元数量
                    act.fct = "logistic", ## 激活函数
                    linear.output = FALSE,
                    algorithm = "rprop+",err.fct ="ce")
###train data
pre_ran<-compute(mlpcla, x_train[,1:23])$net.result
pre_ran <-c(0,1)[apply(pre_ran,1,which.max)]
#net.predict<-ifelse(net.predict == 1,0,1)
#predict.table<-table(y_test,net.predict)
ran_roc <- roc(as.numeric(x_train$class), 
               as.numeric(pre_ran))
###test data
pre_ran2<-compute(mlpcla, x_test[,1:23])$net.result
pre_ran2 <-c(0,1)[apply(pre_ran2,1,which.max)]
ran_roc2 <- roc(as.numeric(x_test$class), 
                as.numeric(pre_ran2))
###SVM model by e1071
###load package
library(caret)
library(e1071)
library(dplyr)
library(pROC)
###split train and test data
x_train<-x_train[,match(ae_marker, colnames(x_train))]
x_test<-x_test[,match(ae_marker, colnames(x_test))]
###SVM model
linear.tune<-tune.svm(class~.,data = x_train,
                      kernel = "linear",
                      cost = c(0.001,0.01,0.1,1,5,10))
best.linear<-linear.tune$best.model
###train data
pre_ran <- predict(best.linear, newdata = x_train[,1:23])
ran_roc <- roc(as.numeric(x_train$class), 
               as.numeric(pre_ran))
###test data
pre_ran2 <- predict(best.linear, newdata = x_test[,1:23])
ran_roc2 <- roc(as.numeric(x_test$class), 
                as.numeric(pre_ran2))
###
###logistic model
model<-glm(class~., data=x_train, family = binomial)
###train data
pre_ran <- predict(model, newdata = x_train[,1:23],type = "response")
pre_ran <- ifelse(pre_ran >= 0.5,1,0)
ran_roc <- roc(as.numeric(x_train$class), 
               as.numeric(pre_ran))
###test data
pre_ran2 <- predict(model, newdata = x_test[,1:23],type = "response")
pre_ran2 <- ifelse(pre_ran2 >= 0.5,1,0)
ran_roc2 <- roc(as.numeric(x_test$class), 
                as.numeric(pre_ran2))
ran_roc2
###XGBoost model
library(xgboost)
library(caret)
library(dplyr)
###XGBoost model
xgb_model <- xgboost(data = as.matrix(x_train[,1:23]), label = as.matrix(x_train[,24]),nrounds = 200, params = params,verbose =0)
###train data
pre_ran <- predict(xgb_model, newdata = as.matrix(x_train[,1:23]))
ran_roc <- multiclass.roc(as.numeric(x_train$class), 
                          as.numeric(pre_ran))
ran_roc
###test data
pre_ran2 <- predict(xgb_model , newdata = as.matrix(x_test[,1:23]))
ran_roc2 <- multiclass.roc(as.numeric(x_test$class), 
                           as.numeric(pre_ran2))
ran_roc2
###randomForest model
library(randomForest)
library(dplyr)
library(pROC)
###
x_train<-x_train[,match(ae_marker, colnames(x_train))]
x_test<-x_test[,match(ae_marker, colnames(x_test))]
###randomForest
set.seed(20230802)
system.time(fit.forest <- randomForest(class~., data=x_train, ntree=3000))
###train data
pre_ran <- predict(fit.forest_2, newdata = x_train)
ran_roc <- roc(as.numeric(x_train$class), 
                          as.numeric(pre_ran))
ran_roc
###test data
pre_ran2 <- predict(fit.forest_2, newdata = x_test)
ran_roc2 <- roc(as.numeric(x_test$class), 
                           as.numeric(pre_ran2))
ran_roc2
