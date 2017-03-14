# Author: Syeda Sakira Hassan
# Date: 17.09.2015
# Leave one out cross validation via glmnet instead of cv.glmnet

leaveOneOutCV <- function(X,y, lambda){
  
  
  # number of samples
  N <- nrow(X)
  
  # initialize
  #lambda <- seq(2^(-3),2^4, by = 0.1)
  model <- vector("list",N)
  
  
  
  corr <- vector()
  mse <- vector()
  mae <- vector()
  
  for( n in 1:N){
    
    # Leave one sample out
    trainIdx <- setdiff(1:N,n)
    testIdx <- n
    
    ## Separate the training and test sets
    Xtrain <- X[trainIdx,]
    Ytrain <- y[trainIdx]
    
    Xtest <- t(X[c(testIdx),])
    Ytest <- y[testIdx]
    
    
    ## train the model with glmnet
    
    fit <- glmnet(Xtrain, Ytrain, lambda = lambda)
    yHat <- predict(fit, Xtest)
    
    model[[n]] <- fit
    #corr[n] <- cor(Ytest, yHat)
    mse[n] <- mean((Ytest-yHat)^2)
    mae[n] <- mean(abs(Ytest-yHat))

  }
  
  results <- list(model = model,  mse = mse, mae = mae)
  
  return(results)
  
  
#   for n = 1:N
#   
#   %% Leave one sample out    
#   trainIdx = setdiff(1:N, n); % all exept one idx
#   testIdx = n;                   % the left out idx
#   
#   Xtrain = X(trainIdx, :);
#   Ytrain = y(trainIdx);
#   
#   
#   
#   
#   P = Xtrain';
#     T = Ytrain';
#   
  
}