# Author: Syeda Sakira Hassan
# Date: 23.07.2014
# ChangeLog Dates: 24.07.2015, 28.07.2015, 30.07.2015, 17.09.2015, 25.11.2015
# 28.12.2015, 31.12.2015, 4.1.2016, 16.1.2017, 17.1.2017, 18.1.2017
# Description: 
# Find the number of genes in the listed pathways for specific bacterias

# clear variables
rm(list = ls())

# load the library
library(KEGGREST)
# the documentation can be found in : 
# http://www.bioconductor.org/packages/release/bioc/vignettes/KEGGREST/inst/doc/KEGGREST-vignette.html


# KEGG bacterial organism code
# updating the orgCodeList by reading the file cdiGMPdataSet2.csv
cdiGMPdataFile <- "P://Users personal data/phd materials/datas/kegg/domains/cdiGMPdataSet2.csv"
cdiGMPdata <- read.csv(cdiGMPdataFile, header=TRUE, sep=';', quote="")
orgCodeList <- cdiGMPdata$ORG.CODE

## Check duplicate values
# n_occur <- data.frame(table(orgCodeList))
# n_occur[n_occur$Freq > 1,]
# orgCodeList[orgCodeList %in% n_occur$Var1[n_occur$Freq > 1],]

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
pathwayList <- sub("path:map", "", refPathways$pathcode)

######## Count the number of genes for each organism in each pathway

# create a matrix
rowNames <- orgCodeList 
colNames <- refPathways$pathcode
orgPathwayCount <- matrix(data = NA, nrow = length(orgCodeList), ncol = length(pathwayList), dimnames = list( rowNames, colNames))


orgPathwaysGenesFilePath <- "P://Users personal data/phd materials/datas/kegg/selected_org_pathways_genes/"
genesPathwayFilePath <-  "P://Users personal data/phd materials/datas/kegg/genes_pathways/"
geneFilePath <- "P://Users personal data/phd materials/datas/kegg/all_gene_list/"





for( i in 1:length(orgCodeList)) {
  
  
  
  orgCode <- orgCodeList[i]
  
  ## read the each organism file and list all the genes
  
  geneFileName <- paste(geneFilePath, orgCode, '.txt', sep="") 
  genes <- read.csv(geneFileName, header=FALSE, sep='\t', quote="")
  colnames(genes) <- c("GENE_ID", "DESCRIPTION")
  
  # create a matrix of all genes and pathways
  rowNames <- genes$GENE_ID
  colNames <- refPathways$pathcode
  
  genesPathwayData <- matrix(data = NA, nrow = nrow(genes), ncol = length(pathwayList), dimnames = list( rowNames, colNames))
  
  
  # read the genes_pathways files
  genesPathwayFile <- paste(genesPathwayFilePath, orgCode, '.txt', sep="")
  genesPathway <- tryCatch({
    read.csv(genesPathwayFile, header=FALSE, sep='\t', quote="")
    
    
    
  },
  error = function(err) {
    
    # error handler picks up where error was generated
    # print(paste("MY_ERROR:  ",err))
    print( paste(err, "\nThe pathways for the organism ", orgCode, " was not found in KEGG database."))
    return(data.frame(matrix(data=NA, ncol = 2)))
  }
  
  
  )# End of tryCatch
  
  
  
  
  colnames(genesPathway) <- c("GENE_ID", "PATHWAY_ID")
  #head(genesPathway)
  
  # search genes in the pathway
  for( j in 1:length(pathwayList)){
    
    searchPathway <- paste("path:", orgCode, pathwayList[j], sep = "")
    selectedGenes <- genesPathway[grep(searchPathway, genesPathway$PATHWAY_ID),]
    genesPathwayData[as.character(selectedGenes$GENE_ID),j] <- 1
    #print(nrow(selectedGenes))
    
    
  }
  
  
  # write to file
  fileName <- paste(orgPathwaysGenesFilePath , orgCode, '.txt', sep="")
  write.table(genesPathwayData,fileName, sep="\t")
  
  totalGenesInPathways <- colSums(genesPathwayData, na.rm = TRUE)
  orgPathwayCount[i,] <- totalGenesInPathways
  
  
  # sum(genesPathwayData[,463], na.rm = TRUE)
  
  
  
}



save(orgPathwayCount, cdiGMPdata, pathwayList, refPathways, file = 'parameters.Rdata')

### Mapping between two different (input and output) files



##################################################################################

# installation guide can be found in here:
# https://cran.r-project.org/web/packages/glmnet/vignettes/glmnet_beta.html
#install.packages("glmnet", repos = "http://cran.us.r-project.org")
library(glmnet)

x <- orgPathwayCount
yPilZdomain <- cdiGMPdata$PilZ

#Normalized Data
normalized = (x-min(x))/(max(x)-min(x))

#Histogram of example data and normalized data
# par(mfrow=c(1,2))
# hist(x,xlab="Data",col="lightblue",main="")
# hist(normalized,xlab="Normalized Data",col="lightblue",main="")



n <- 5 # 5-fold cross validation

cvfit = cv.glmnet(x, yPilZdomain, nfolds = n, standardize.response = TRUE)


plot(cvfit)
coef(cvfit, s = "lambda.min")
cvfit


###################################
###################################
### start from here glmnet
fit <- glmnet(x,yPilZdomain)
fit
plot(fit)
predict(fit,newx=x[1:5,],s=c(0.01,0.005))

# defualt lambdas
lambdas <- fit$lambda
source("P:/Users personal data/phd materials/datas/kegg/getBinaryBeeError.R")
X <- x
y <- yPilZdomain

errBEE <- 0

for(l in 1:length(lambdas)){
  lamb <- lambdas[l]
  a <- fit$beta[,l]
  b <- fit$a0[l]
  
  err <- getBinaryBeeError(X,y,a,b)
  errBEE[l] <- err$err
  print(paste0( "For lambda = ", lamb, " BEE error: ", err$err))
}

plot(lambdas, errBEE, type='l')
save(X,y, lambdas, fit, errBEE, file = 'modelParams.Rdata')

set.seed(1010)
cvfit = cv.glmnet(x, y,   nfolds = 5)
#cvfitNorm = cv.glmnet(normalized, y, lambda = lambda, nfolds = 3)
cvfit
#cvfitNorm
plot(cvfit)
yHat1 <- predict(cvfit, newx = x, s = 'lambda.min')
yHat2 <- predict(cvfit, newx = x, s = 'lambda.1se')
cor( y, yHat1)
cor(y, yHat2)
yHat1
yHat2

set.seed(1010)
lambda = seq(2^(-3),2^4, by = 0.1)
cvfit = cv.glmnet(x,y, lambda = lambda, nfolds = 3)
cvfit
plot(cvfit)
yHat3 <- predict(cvfit,newx=x, s="lambda.min")
yHat4 <- predict(cvfit,newx=x, s="lambda.1se")

cor( y, yHat1) # correlation 0.89%
cor( y, yHat2)
cor( y, yHat3)
cor( y, yHat4)

########## PilZdomain

set.seed(101)
cvfit = cv.glmnet(x, yPilZdomain,   nfolds = 5)
cvfit
plot(cvfit)
yHat1 <- predict(cvfit, newx = x, s = 'lambda.min')
cor( yPilZdomain, yHat1)
yHat1


cor( yPilZdomain, yHat1) # correlation 0.89%


coef <- coef(cvfit, s = "lambda.min")
coef
beta  <- as.vector( t(coef( cvfit ,s="lambda.min"))) 
beta
testX <- cbind(1,x) 
yhat2  <- testX %*% beta 
plot(yPilZdomain,yhat2)
plot(cvfit)
cbind(yPilZdomain,yhat2, (yPilZdomain-yhat2)^2)
sum((yPilZdomain-yhat2)^2)

coefList <- data.frame("pathcode" = "Intercept", "description" = "Intercept")
coefList <- rbind(coefList,refPathways)
coefList <- cbind(coefList, V3 = NA)
coefList$V3 <- beta
coefList

#################### 
########## Random Forest
#install.packages('randomForest')
library(randomForest)
# Documentation Links:
# http://trevorstephens.com/post/73770963794/titanic-getting-started-with-r-part-5-random
# https://cran.r-project.org/web/packages/randomForest/randomForest.pdf
set.seed(10)
fit <- randomForest(x,yPilZdomain, ntree = 100)
fit
varImpPlot(fit)
pred <- predict(fit, x)
cbind(yPilZdomain, pred)
cor(yPilZdomain, pred)
importance(fit)

##########################
#### PLS
# link: http://www.r-bloggers.com/partial-least-squares-regression-in-r/
install.packages('pls')
library(pls)
#data <- cbind(x,y)
#fitPLS <- plsr(x,yPilZdomain)
install.packages('plsdepot')
library(plsdepot)
fitPLS <- plsreg1(x, yPilZdomain, comps = 3)
fitPLS
plot(fitPLS)
#plot each observation predicted versus actual
plot(yPilZdomain, fitPLS$y.pred, type = "n", xlab="Original", ylab = "Predicted")
title("Comparison of responses", cex.main = 0.9)
abline(a = 0, b = 1, col = "gray85", lwd = 2)
text(yPilZdomain, fitPLS$y.pred, col = "#5592e3")
###################################################################
####################################################################
############ GLMNET LOOP

########## PilZdomain
corr <- vector()
mse <- vector()
for (i in 1:10){
  set.seed(i)
  cvfit = cv.glmnet(x, yPilZdomain,   nfolds = 5)
  plot(cvfit)
  yHat1 <- predict(cvfit, newx = x, s = 'lambda.min')
  corr[i] <- cor( yPilZdomain, yHat1)
  
  
  beta  <- as.vector( t(coef( cvfit ,s="lambda.min"))) 
  testX <- cbind(1,x) 
  yhat2  <- testX %*% beta 
  plot(yPilZdomain,yhat2, type = 'l')
  
  mse[i] <- sum((yPilZdomain-yhat2)^2)
  
  coefList <- data.frame("pathcode" = "Intercept", "description" = "Intercept")
  coefList <- rbind(coefList,refPathways)
  coefList <- cbind(coefList, V3 = NA)
  coefList$V3 <- beta
  cat('Iteration ',i)
  cat('Correlation = ', corr[i], ' MSE = ', mse[i])
  print(coefList)
}
cbind(corr, mse)
numberOfOrg <- length(yPilZdomain)

plot(1:numberOfOrg,yPilZdomain, type = 'l', cex=0.1, col = "blue", xlab = 'Representative complete genomes from all bacterial species', ylab = 'Number of genes encoding PilZ domain')
points(1:numberOfOrg,yPilZdomain, cex = .5, col = "dark blue")

lines(1:numberOfOrg,yhat2, col = 'red', cex=0.1)
points(1:numberOfOrg,yhat2, cex = .5, col = "dark red")
legend(5, 30, c("True", "Predicted"), col = c('blue', 'red'),
       text.col = "green4", lty = c(1, 1), pch = c(19, 19),
       merge = TRUE, bg = "gray90") # lty = line type, pch = point type on the line
#####################################################
###################################################
#### Random Forest Loop

corrRF <- vector()
mseRF <- vector()
for (i in 1:10){
  
  set.seed(i)
  fit <- randomForest(x,yPilZdomain, ntree = 100)
  fit
  varImpPlot(fit)
  pred <- predict(fit, x)
  
  corrRF[i] <- cor(yPilZdomain, pred)
  mseRF[i] <- sum((yPilZdomain-pred)^2)
  impRF <- importance(fit)
  
  coefListRF <- data.frame("pathcode" = character(), "description" = character())
  coefListRF <- rbind(coefListRF,refPathways)
  coefListRF <- cbind(coefListRF, V3 = NA)
  coefListRF$V3 <- impRF
  
  cat('Iteration ',i)
  cat('Correlation = ', corrRF[i], ' MSE = ', mseRF[i])
  print(coefListRF)
  
  
}
cbind(corrRF, mseRF)

plot(1:numberOfOrg,yPilZdomain, type = 'l', cex = .1, col = "blue", xlab = 'Representative complete genomes from all bacterial species', ylab = 'Number of genes encoding PilZ domain')
points(1:numberOfOrg,yPilZdomain, cex = .5, col = "dark blue")

lines(1:numberOfOrg,pred, col = 'red', cex = .1)
points(1:numberOfOrg,pred, cex = .5, col = "dark red")
legend(5, 30, c("True", "Predicted"), col = c('blue', 'red'),
       text.col = "green4", lty = c(1, 1), pch = c(19, 19),
       merge = TRUE, bg = "gray90") # lty = line type, pch = point type on the line

#############################################################
############ COMPARISON PLOTS
plot(1:numberOfOrg,yPilZdomain, type = 'l', col = "blue", xlab = 'indices', ylab = 'prediction')
points(1:numberOfOrg,yPilZdomain, cex = .5, col = "dark blue")

lines(1:numberOfOrg,pred, col = 'red')
points(1:numberOfOrg,pred, cex = .5, col = "dark red")
lines(1:numberOfOrg,yhat2, col = 'green')
points(1:numberOfOrg,yhat2, cex = .5, col = "dark green")
legend(8.5, 5, c("True", "prediction RF (mse = 5.03, corr = 98%)", "prediction GLMnet (mse = 3.19, corr = 95%)"), col = c('blue', 'red', 'green'),
       text.col = "green4", cex=0.5, lty = c(1, 1, 1), pch = c(19, 19, 19),
       merge = TRUE, bg = "gray90") # lty = line type, pch = point type on the line


# Example of adding function script and obtain the result
# setwd("P:/Users personal data/phd materials/datas/kegg")
# source("addPercent.R")
# x <- c(0.458, 1.6653, 0.83112)
# print(addPercent(x))


##############################################################
## Leave One Out CV
# source('loo_cross_val.R')
# res <- leaveOneOutCV(x,yPilZdomain)
# res
# res$mse
# 
# (res$mse)/(length(res$mse))

