# script for measuring Mean squared error
# With 5-fold cross-validation 
# With leave-one-out cross-validation
# With Bayesian Error estimation
# ChangeLog: 04.01.2016, 11.01.2016, 23.1.2017, 24.1.2017

rm(list = ls())
# set the working directory
dir <- "P:/Users personal data/phd materials/datas/kegg"
setwd(dir = dir)
getwd()

# load the data
load("parameters.Rdata")

# load the library
# for glmnet
library(glmnet)
# for random forest
library(randomForest)


pathwayNames <- colnames(orgPathwayCount)
orgNames <- rownames(orgPathwayCount)


# total number of samples
N <- nrow(orgPathwayCount)
numberOfDomains <- 5

X <- orgPathwayCount
# y is the matrix of c-di-GMP signaling domains gene counts
y <- matrix(data = NA, nrow = N, ncol = numberOfDomains)
y[,1] <- cdiGMPdata$GGDEF.only
y[,2] <- cdiGMPdata$EAL.only
y[,3] <- cdiGMPdata$HD.GYP
y[,4] <- cdiGMPdata$PilZ
y[,5] <- cdiGMPdata$MshEN

colnames(y) <- c("GGDEF", "EAL", "HD-GYP", "PilZ", "MshEN")






# initialize output variables
#yHatGLMNET <- vector(length = 5)
yHatGLMNET <- matrix(data = NA, nrow = N, ncol = numberOfDomains)
yHatRF <- matrix(data = NA, nrow = N, ncol = numberOfDomains)


n <- 1
for(n in 1:N){
  
  # Leave one sample out
  trainIdx <- setdiff(1:N,n)
  testIdx <- n
  
  
  orgNameTest <- orgNames[n]
  
  
  for (i in 1:numberOfDomains){
    ## Separate the training and test sets
    Xtrain <- X[trainIdx,]
    Ytrain <- y[trainIdx,i]
  
    Xtest <- t(X[c(testIdx),])
    Ytest <- y[testIdx,i]
  
  
  
    #   Xtest
    #   Ytest
    #   orgNameTest
  
  
    # 10-fold cross validation glmnet model
  
    # number of folds
    nfolds <- 10 
    modelGLMNET <- cv.glmnet(Xtrain, Ytrain, nfolds = nfolds, standardize.response = TRUE)
    yHatGLMNET[n,i] <- predict(modelGLMNET, newx = Xtest, s = "lambda.min")
    # yHatGLMNET
  
    # Random forest model
    set.seed(n)
    modelRF <- randomForest(Xtrain, Ytrain, proximity = TRUE)
    yHatRF[n,i] <- predict(modelRF, Xtest, type = "response")
    
  }
  
  print(cat("Prediction for", orgNameTest , "organism is completed...\n"))
  
  

}

getwd()
save(yHatGLMNET, yHatRF, file = 'results.Rdata')

# Plotting

# plot(1:N, y, type = 'l', cex = .1, col = "blue", xlab = 'Representative complete genomes from all bacterial species', ylab = 'Number of genes encoding PilZ domain')
# points(1:N, y, cex = .5, col = "dark blue")
# 
# lines(1:N, yHatGLMNET, col = 'red', cex = .1)
# points(1:N, yHatGLMNET, cex = .5, col = "dark red")
# legend(5, 30, c("True", "Predicted"), col = c('blue', 'red'),
#        text.col = "green4", lty = c(1, 1), pch = c(19, 19),
#        merge = TRUE, bg = "gray90") # lty = line type, pch = point type on the line
# 
# plot(1:N, y, type = 'l', cex = .1, col = "blue", xlab = 'Representative complete genomes from all bacterial species', ylab = 'Number of genes encoding PilZ domain')
# points(1:N, y, cex = .5, col = "dark blue")
# 
# lines(1:N, yHatRF, col = 'red', cex = .1)
# points(1:N, yHatRF, cex = .5, col = "dark red")
# legend(5, 30, c("True", "Predicted"), col = c('blue', 'red'),
#        text.col = "green4", lty = c(1, 1), pch = c(19, 19),
#        merge = TRUE, bg = "gray90") # lty = line type, pch = point type on the line
# 
# # sorting
# Ysorted <- sort(y, index.return = TRUE)
# plot(1:N, Ysorted$x, type = 'p', cex = .1, col = "blue", xlab = 'Representative complete genomes from all bacterial species', ylab = 'Number of genes encoding PilZ domain')
# points(1:N, Ysorted$x, cex = .5, col = "dark blue")
# 
# yHatRFsorted <- yHatRF[Ysorted$ix]
# lines(1:N, yHatRFsorted, col = 'red', cex = .1)
# points(1:N, yHatRFsorted, cex = .5, col = "dark red")
# legend(5, 30, c("True", "Predicted"), col = c('blue', 'red'),
#        text.col = "green4", lty = c(1, 1), pch = c(19, 19),
#        merge = TRUE, bg = "gray90") # lty = line type, pch = point type on the line
# 

# setEPS()
# postscript("whatever.eps")
# plot(rnorm(100), main="Hey Some Data")
# dev.off()




domainNames <- colnames(y)

# Plotting and save as pdf
i <- 1
for(i in 1:numberOfDomains){
  
  
  # sorting
  Ysorted <- sort(y[,i], index.return = TRUE)
  yHatGLMNETsorted <- yHatGLMNET[Ysorted$ix, i]
  yHatRFsorted <- yHatRF[Ysorted$ix, i]
  
  
  # for Lasso
  fileName <- paste("images/Lasso", domainNames[i], ".pdf", sep="" )
  pdf(fileName)
  plot(1:N, Ysorted$x, type = 'p', cex = .1,  col = "blue", xlab = 'Representative complete genomes from all bacterial species', ylab = paste('Number of genes encoding ', domainNames[i], ' domain', sep="") )
  
  points(1:N, yHatGLMNETsorted, cex = .5, pch = 4, col = "dark red")
  legend(5, 30, c("True", "Predicted"), col = c('blue', 'red'),
         text.col = "green4", lty = c(1, 1), pch = c(19, 4),
         merge = TRUE, bg = "gray90") # lty = line type, pch = point type on the line
  
  
  
  dev.off()
  
  # for Random Forest
  fileName <- paste("images/RandomForest", domainNames[i], ".pdf", sep="" )
  pdf(fileName)
  plot(1:N, Ysorted$x, type = 'p', cex = .1, col = "blue", xlab = 'Representative complete genomes from all bacterial species', ylab = paste('Number of genes encoding ', domainNames[i], ' domain', sep="") )
  
  points(1:N, yHatRFsorted, cex = .5, pch = 4, col = "dark red")
  legend(5, 30, c("True", "Predicted"), col = c('blue', 'red'),
         text.col = "green4", lty = c(1, 1), pch = c(19, 4),
         merge = TRUE, bg = "gray90") # lty = line type, pch = point type on the line
  
  dev.off()
}


## Finding significant pathways 

numberOfPathways <- 12
coefGLMNET <- matrix(data = NA, nrow = numberOfPathways, ncol = numberOfDomains)
coefRF <- matrix(data = NA, nrow = numberOfPathways, ncol = numberOfDomains)

i <- 1
for ( i in 1:numberOfDomains){
  
  # modeling
  
  set.seed(100)
  modelGLMNET <- cv.glmnet(X, y[,i], nfolds = nfolds, standardize.response = TRUE)
  
  modelRF <- randomForest(X, y[,i], proximity = TRUE)
  
  print(cat ("modeling for ", domainNames[i], " is done..."))
  
  # save coefficients
  coefficients <- coef(modelGLMNET)
  coefGLMNET[,i] <- coefficients[2:length(coefficients)]
  
  coefficients <- importance(modelRF)
  coefficients <- scale(coefficients, center = TRUE, scale = TRUE)
  coefRF[,i] <- coefficients
}
  
print("Done...")

# plot heatmap

# prepare the data
coefGLMNET[coefGLMNET==0] <- NA

nba_matrix <- data.matrix(cbind(coefGLMNET, coefRF))

labCol <- paste("Lasso(", domainNames, ")", sep="")
labCol <- c(labCol, paste("Random Forest(", domainNames, ")", sep=""))

labRow <- data.frame(pathcode = c("path:map00010", "path:map00020", "path:map00030", "path:map00500", "path:map00520",
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


# make a heatmap
#library(gplots)
#install.packages("pheatmap")
# install.packages("RColorBrewer")

library(RColorBrewer)

library(pheatmap)
pheatmap(nba_matrix, 
         col = greenred(256), 
         cluster_cols = F, 
         fontsize_row = 2, 
         border_color = NA)
pdf("images/PathwaysHeatmap.pdf")
nba_heatmap <- heatmap(nba_matrix, 
                       Rowv=NA, 
                       Colv = NA,  
                       labRow = labRow$description, 
                       labCol = labCol,  
                       col = brewer.pal(9, "PuBu"),   
                       scale="column", 
                       cexRow = 0.6,
                       cexCol = 0.8,
                       
                       margins=c(10,10))
dev.off()