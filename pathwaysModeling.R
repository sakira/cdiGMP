#save(celluloseCountGenome, file = 'celluloseCountGenome.Rda')
#save(genomesMetaInfo, file = 'genomesMetaInfo.Rda')
rm(list = ls())

pathname <- "P://Users personal data/phd materials/datas/imgdb/bacterialcellulose/"
setwd(pathname)

load('genomesMetaInfo2.Rda')
load('celluloseCountGenome2.Rda')

# for sorting
sortedCelluloseCountGenome <- celluloseCountGenome[order(celluloseCountGenome$Genome, decreasing = TRUE),]

# combine the number of genes in pathways and cellulose 
fileName <- "P://Users personal data/phd materials/datas/kegg/geneCountInPathways.txt"
genesInPathways <- read.csv(fileName, header=TRUE, sep='\t')

# remove duplicate just sake of easy computation
genesInPathways <- genesInPathways[!duplicated(genesInPathways$organism),]
sortedgenesInPathways <- genesInPathways[order(genesInPathways$organism, decreasing = TRUE),]
xSum <- rowSums(as.matrix(sortedgenesInPathways[c(3:465)]))
sortedgenesInPathways <- cbind(sortedgenesInPathways, genes = 0)
sortedgenesInPathways$genes[1:2731] <- xSum

# add new column

genesInPathways <- cbind(genesInPathways, cellulose = 0)
#cellulose <- data.frame(matrix(nrow=1, ncol = ))

i <- 1

for (org in genesInPathways){
  
  tempOrg <- genesInPathways$organism[i]
  
  tempInfo <- celluloseCountGenome[celluloseCountGenome$Genome == tempOrg,]
  
  
  if(nrow(tempInfo) != 0 ){
    genesInPathways$cellulose[i] <- tempInfo$Gene.Count
  }
  i <- i + 1
  
}




x <- as.matrix(genesInPathways[c(3:465)])
y <- genesInPathways[c(466)]
lmfit <- lm(x  ~  y$cellulose)
sink('outputLM2.txt')
summary(lmfit)
sink()
glmfit <- glm(y$cellulose ~ x)
sink('outputGLM2.txt')
summary(glmfit)
sink()

barplot(coef(glmfit))

coefficients <- coef(glmfit)


## read the reference pathway file & list them
fileName <- "P://Users personal data/phd materials/datas/kegg/ref_pathways.txt"
refPathways <- read.csv(fileName, header=FALSE, sep='\t')

coefOfPathways <- data.frame(V1 = "Intercept", V2 = "Intercept")
coefOfPathways <- rbind(coefOfPathways, refPathways)
coefOfPathways <- cbind(coefOfPathways, V3 = NA)
coefOfPathways$V3 <- coefficients
barplot(coefOfPathways$V3, names.arg = coefOfPathways$V2)

X <- summary.glm(glmfit)

sink('ListOfPathwaysCoef2.txt')
coefOfPathways
sink()

#X <- coef(summary(glmfit))
coefOfPathways <- cbind(coefOfPathways, V4 = NA)
#coefOfPathways$V4[1:nrow(X)] <- X$coefficients[,4] # get the probability values
coefSignificance <- X$coefficients[,4] # get the probability values
colNames <- coefOfPathways$V1[which(!is.na(coefOfPathways$V3))]

#coefOfPathways$V4[colNames] <- coefSignificance[colNames]
col <- (!is.na(coefOfPathways$V3))
coefOfPathways$V4[col] <- coefSignificance
sink('ListOfPathwaysSignificance2.txt')
coefOfPathways
sink()



X <- order(coefOfPathways$V4)
Y <- coefOfPathways[X,]
View(Y)
sink('SortedListOfPathwaysSignificance2.txt')
Y
sink()


colnames(Y) <- c("KEGG Pathway ID", "Pathway description", "Coefficients", "Significance Pr(>|t|)")



require(graphics)
require(grDevices)
x  <- as.matrix(genesInPathways[c(3:466)])
x <- as.matrix(genesInPathways[c(1:1000),c(3:100,466)])
rc <- rainbow(nrow(x), start = 0, end = .3)
cc <- rainbow(ncol(x), start = 0, end = .3)

org <- genesInPathways[c(1:1000),c(2)]


hv <- heatmap(x, Colv = NA, Rowv = NA, col = cm.colors(256), scale = "column",
              RowSideColors = rc, ColSideColors = cc, 
              labCol = refPathways$V2[c(1:99)], labRow = org,
              main = "heatmap" )
#utils::str(hv) # the two re-ordering index vectors



library(gplots)
heatmap.2(x, Rowv=FALSE, Colv=FALSE, 
          dendrogram="none", col = c("blue", "red"), RowSideColors = rc, ColSideColors = cc, 
          key=T, keysize=1.5, density.info="none", 
          trace="none", labCol = refPathways$V2[c(1:99)], labRow=org)














