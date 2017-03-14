# script for measuring Mean squared error
# With 5-fold cross-validation 
# With leave-one-out cross-validation
# With Bayesian Error estimation
# ChangeLog: 04.01.2016, 11.01.2016

rm(list = ls())
# set the working directory
dir <- "P:/Users personal data/phd materials/datas/kegg"
setwd(dir = dir)
getwd()

# load the data
load("modelParams.Rdata")


# bayesian error estimation
source("P:/Users personal data/SOFTWARES/BEE/getBinaryBeeError.R")

#fit <- glmnet(X, y, family = "binomial")
#lambdas <- fit$lambda

errBEE <- 0


for(l in 1:length(lambdas)){
  lamb <- lambdas[l]
  a <- fit$beta[,l]
  b <- fit$a0[l]
  
  err <- getBinaryBeeError(X,y,a,b)
  errBEE[l] <- err$err
  print(paste0( "For lambda = ", lamb, " BEE error: ", err$err))
}


plot(log10(lambdas), errBEE)
lines(log10(lambdas),errBEE, type = 'l')
# Find the lambda which gives minimum error
minVal <- min(errBEE)
lambdaIdx <- which.min(errBEE)
minLambda <- lambdas[lambdaIdx]

library(glmnet)

# prediction for the minimum lambda
yHat <- predict(fit, newx = X, s = minLambda, type = "response")

corrBEE <- cor( y, yHat)
mseBEE <- mean((y-yHat)^2)
maeBEE <- mean(abs(y-yHat))
corrBEE
mseBEE
maeBEE

###############################################################
#### 5-fold cross-validation approach
corr <- vector()
mse <- vector()
mae <- vector()

for (i in 1:100){
  set.seed(i)
  cvfit = cv.glmnet(X, y, lambda = lambdas,  nfolds = 5)
  yHat <- predict(cvfit, newx = X, s = 'lambda.min')
  
  corr[i] <- cor( y, yHat)
  mse[i] <- mean((y-yHat)^2)
  mae[i] <- mean(abs(y-yHat))

  
  cat('Iteration ',i)
  cat(' Correlation = ', corr[i], ' MSE = ', mse[i], '\n')
  
}

mseCV5 <- mean(mse)
corrCV5 <- mean(corr)
maeCV5 <- mean(mae)

corrCV5
mseCV5
maeCV5


###############################################################
#### leave-one-out cross-validation approach

source('loo_cross_val.R')
res <- leaveOneOutCV(X,y, lambdas)
mseLOO <- mean(res$mse)
maeLOO <- mean(res$mae)
#corrLOO <- res$corr

mseLOO
maeLOO


#########################################################
#########################################################
# Other domains for output
cdiGMPdataFile <- "P://Users personal data/phd materials/datas/kegg/domains/cdiGMPdataSet1.csv"
cdiGMPdata <- read.csv(cdiGMPdataFile, header=TRUE, sep=';', quote="")

y <- cdiGMPdata$GGDEF.only

# bayesian error estimation
source("P:/Users personal data/SOFTWARES/BEE/getBinaryBeeError.R")

fit <- glmnet(X, y, family = "binomial")
lambdas <- fit$lambda

errBEE <- 0


for(l in 1:length(lambdas)){
  lamb <- lambdas[l]
  a <- fit$beta[,l]
  b <- fit$a0[l]
  
  err <- getBinaryBeeError(X,y,a,b)
  errBEE[l] <- err$err
  print(paste0( "For lambda = ", lamb, " BEE error: ", err$err))
}


plot(log10(lambdas), errBEE)
lines(log10(lambdas),errBEE, type = 'l')

# Find the lambda which gives minimum error
minVal <- min(errBEE)
lambdaIdx <- which.min(errBEE)
minLambda <- lambdas[lambdaIdx]

library(glmnet)

# prediction for the minimum lambda
yHat <- predict(fit, newx = X, s = minLambda, type = "response")

corrBEE <- cor( y, yHat)
mseBEE <- mean((y-yHat)^2)
maeBEE <- mean(abs(y-yHat))

corrBEE
mseBEE
maeBEE

###############################################################
#### 5-fold cross-validation approach
corr <- vector()
mse <- vector()
mae <- vector()
for (i in 1:100){
  set.seed(i)
  cvfit = cv.glmnet(X, y, lambda = lambdas,  nfolds = 5)
  yHat <- predict(cvfit, newx = X, s = 'lambda.min')
  
  corr[i] <- cor( y, yHat)
  mse[i] <- mean((y-yHat)^2)
  mae[i] <- mean(abs(y-yHat))

  
  cat('Iteration ',i)
  cat(' Correlation = ', corr[i], ' MSE = ', mse[i], '\n')
  
}

mseCV5 <- mean(mse)
corrCV5 <- mean(corr)
maeCV5 <- mean(mae)
corrCV5
mseCV5
maeCV5



###############################################################
#### leave-one-out cross-validation approach

source('loo_cross_val.R')
res <- leaveOneOutCV(X,y, lambdas)
mseLOO <- mean(res$mse)
maeLOO <- mean(res$mae)
#corrLOO <- res$corr

mseLOO
maeLOO
##########################################
###########################################
########################################
#########################################

y <- cdiGMPdata$GGDEF.EAL

# bayesian error estimation
source("P:/Users personal data/SOFTWARES/BEE/getBinaryBeeError.R")

fit <- glmnet(X, y, family = "binomial")
lambdas <- fit$lambda

errBEE <- 0


for(l in 1:length(lambdas)){
  lamb <- lambdas[l]
  a <- fit$beta[,l]
  b <- fit$a0[l]
  
  err <- getBinaryBeeError(X,y,a,b)
  errBEE[l] <- err$err
  print(paste0( "For lambda = ", lamb, " BEE error: ", err$err))
}


plot(log10(lambdas), errBEE)
lines(log10(lambdas),errBEE, type = 'l')

# Find the lambda which gives minimum error
minVal <- min(errBEE)
lambdaIdx <- which.min(errBEE)
minLambda <- lambdas[lambdaIdx]

library(glmnet)

# prediction for the minimum lambda
yHat <- predict(fit, newx = X, s = minLambda, type = "response")

corrBEE <- cor( y, yHat)
mseBEE <- mean((y-yHat)^2)
maeBEE <- mean(abs(y-yHat))

corrBEE
mseBEE
maeBEE
###############################################################
#### 5-fold cross-validation approach
corr <- vector()
mse <- vector()
mae <- vector()
for (i in 1:100){
  set.seed(i)
  cvfit = cv.glmnet(X, y, lambda = lambdas,  nfolds = 5)
  yHat <- predict(cvfit, newx = X, s = 'lambda.min')
  
  corr[i] <- cor( y, yHat)
  mse[i] <- mean((y-yHat)^2)
  mae[i] <- mean(abs(y-yHat))

  
  cat('Iteration ',i)
  cat(' Correlation = ', corr[i], ' MSE = ', mse[i], '\n')
  
}

mseCV5 <- mean(mse)
corrCV5 <- mean(corr)
maeCV5 <- mean(mae)
corrCV5
mseCV5
maeCV5



###############################################################
#### leave-one-out cross-validation approach

source('loo_cross_val.R')
res <- leaveOneOutCV(X,y, lambdas)
mseLOO <- mean(res$mse)
maeLOO <- mean(res$mae)
#corrLOO <- res$corr

mseLOO
maeLOO

##########################################
###########################################
########################################
#########################################


y <- cdiGMPdata$EAL.only

# bayesian error estimation
source("P:/Users personal data/SOFTWARES/BEE/getBinaryBeeError.R")

fit <- glmnet(X, y, family = "binomial")
lambdas <- fit$lambda

errBEE <- 0


for(l in 1:length(lambdas)){
  lamb <- lambdas[l]
  a <- fit$beta[,l]
  b <- fit$a0[l]
  
  err <- getBinaryBeeError(X,y,a,b)
  errBEE[l] <- err$err
  print(paste0( "For lambda = ", lamb, " BEE error: ", err$err))
}


plot(log10(lambdas), errBEE)
lines(log10(lambdas),errBEE, type = 'l')

# Find the lambda which gives minimum error
minVal <- min(errBEE)
lambdaIdx <- which.min(errBEE)
minLambda <- lambdas[lambdaIdx]

library(glmnet)

# prediction for the minimum lambda
yHat <- predict(fit, newx = X, s = minLambda, type = "response")

corrBEE <- cor( y, yHat)
mseBEE <- mean((y-yHat)^2)
corrBEE
mseBEE

###############################################################
#### 5-fold cross-validation approach
corr <- vector()
mse <- vector()
for (i in 1:100){
  set.seed(i)
  cvfit = cv.glmnet(X, y, lambda = lambdas,  nfolds = 5)
  yHat <- predict(cvfit, newx = X, s = 'lambda.min')
  
  corr[i] <- cor( y, yHat)
  mse[i] <- mean((y-yHat)^2)
  
  
  cat('Iteration ',i)
  cat(' Correlation = ', corr[i], ' MSE = ', mse[i], '\n')
  
}

mseCV5 <- mean(mse)
corrCV5 <- mean(corr)
corrCV5
mseCV5


###############################################################
#### leave-one-out cross-validation approach

source('loo_cross_val.R')
res <- leaveOneOutCV(X,y, lambdas)
mseLOO <- mean(res$mse)
#corrLOO <- res$corr

mseLOO

##########################################
###########################################
########################################
#########################################


###################################################
#### Random Forest approach

library(randomForest)
corr <- vector()
mse <- vector()
for (i in 1:100){
  
  set.seed(i)
  fit <- randomForest(X, y, ntree = 100)
  
  
  pred <- predict(fit, X)
  corr[i] <- cor(y, pred)
  mse[i] <- mean((y - pred)^2)
  #impRF <- importance(fit)
  
  
  
  cat('Iteration ',i)
  cat(' Correlation = ', corr[i], ' MSE = ', mse[i], '\n')
  
  
  
}

corrRF <- mean(corr)
mseRF <- mean(mse)
corrRF
mseRF



################################################################
## Reference pathways which we are interested in or we think are
# significant 

refPathways <- data.frame(pathcode = c("path:map00010", "path:map00020", "path:map00030", "path:map00500", "path:map00520",
                                       "path:map00540", "path:map00600", "path:map00900", "path:map01230", "path:map02010",
                                       "path:map02030", "path:map02060"), 
                          description = c( "Glycolysis / Gluconeogenesis",
                                           "Citrate cycle (TCA cycle)",
                                           "Pentose phosphate pathway",
                                           "Starch and sucrose metabolism",
                                           "Amino sugar and nucleotide sugar metabolism",
                                           "Lipopolysaccharide biosynthesis",
                                           "Sphingolipid metabolism",
                                           "Terpenoid backbone biosynthesis",
                                           "Biosynthesis of amino acids",
                                           "ABC transporters",
                                           "Bacterial chemotaxis",
                                           "Phosphotransferase system (PTS)"
                          ))


###################################################################
####################################################################
############ GLMNET LOOP

# library(glmnet)
# ########## PilZdomain
# corr <- vector()
# mse <- vector()
# for (i in 1:10){
#   set.seed(i)
#   cvfit = cv.glmnet(X, y, lambda = lambdas,  nfolds = 5)
#   plot(cvfit)
#   yHat1 <- predict(cvfit, newx = X, s = 'lambda.min')
#   corr[i] <- cor( y, yHat1)
#   
#   
#   beta  <- as.vector( t(coef( cvfit ,s="lambda.min"))) 
#   testX <- cbind(1,X) 
#   yhat2  <- testX %*% beta 
#   plot(y,yhat2, type = 'l')
#   
#   mse[i] <- mean(((y-yhat2)^2))
#   
#   coefList <- data.frame("pathcode" = "Intercept", "description" = "Intercept")
#   coefList <- rbind(coefList,refPathways)
#   coefList <- cbind(coefList, V3 = NA)
#   coefList$V3 <- beta
#   cat('Iteration ',i)
#   cat('Correlation = ', corr[i], ' MSE = ', mse[i])
#   print(coefList)
# }
# 
# cbind(corr, mse)
# 
# n <- length(y)
# 
# plot(1:n,y, type = 'p', col = "blue", xlab = 'indices', ylab = 'prediction')
# points(1:n,y, cex = .5, col = "dark blue")
# 
# lines(1:n,yhat2, col = 'red')
# points(1:n,yhat2, cex = .5, col = "dark red")
# legend(2, 5, c("PilZdomain", "yhat"), col = c('blue', 'red'),
#        text.col = "green4", lty = c(1, 1), pch = c(19, 19),
#        merge = TRUE, bg = "gray90") # lty = line type, pch = point type on the line
