# Author: Syeda Sakira Hassan
# Date: 19.11.2014
# ChangeLog Dates: 08.08.2014

# Description: 


rm(list = ls())



## read the number of genes available in reference pathway files for each bacterial genome/species

fileName <- "P://Users personal data/phd materials/datas/kegg/geneCountInPathways.txt"
genesInPathways <- read.csv(fileName, header=FALSE, sep='\t')
head(genesInPathways)


#fileName <- "P://Users personal data/phd materials/datas/kegg/organism.txt"
#organisms <- read.csv(fileName, header=FALSE, sep='\t')
#head(organisms)

# combine all the cellulose producing genes of 
# different bacterial genomes from all files in the directory
# celluloseOrg variable contails all the listed cellulose genes and genomes



pathname <- "P://Users personal data/phd materials/datas/imgdb/bacterialcellulose/"
setwd(pathname)
fileList <- list.files()
View(fileList)


for(file in fileList){
  # if the merged dataset doesn't exist, create it
  if (!exists("celluloseOrg")){
    celluloseOrg <- read.csv(file, header=TRUE, sep="\t")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("celluloseOrg")){
    temp_dataset <- read.csv(file, header=TRUE, sep="\t")
    celluloseOrg <- rbind(celluloseOrg, temp_dataset)
    rm(temp_dataset)
  }
  
  
}


# genelist19797_06-feb-2015.xls files contains the genes for EC 2.4.1.12
file <- "P://Users personal data/phd materials/datas/imgdb/bacterialcellulose/genelist19797_06-feb-2015.xls"
celluloseOrg <- read.csv(file, header=TRUE, sep="\t")



# find the unique genomes from the celluloseOrg variables
genomes <- unique(celluloseOrg$Genome)

# create a empty data frame
celluloseCountGenome <- data.frame(matrix(ncol=5, nrow=length(genomes)))
colnames(celluloseCountGenome) <- c( "Genome", "Gene.Count", "Gene.ID", "Locus.Tag", "Gene.Product.Name" )

i <- 1

for ( org in genomes){
  #print(org)
  # find the genome in 'celluloseOrg'
  tempOrg <- celluloseOrg[celluloseOrg$Genome.Name == org,]
  #print(tempOrg)
  # save the values in data frame
  celluloseCountGenome[i,]$Genome <- org
  celluloseCountGenome[i,]$Gene.Count <- nrow(tempOrg)
  celluloseCountGenome[i,]$Gene.ID <- paste((tempOrg$Gene.ID),  collapse= ";")
  celluloseCountGenome[i,]$Locus.Tag <- paste((tempOrg$Locus.Tag),  collapse= ";")
  celluloseCountGenome[i,]$Gene.Product.Name <- paste((tempOrg$Gene.Product.Name),  collapse= ";")
  
  rm(tempOrg)
  
  i <- i + 1
  
  
}


celluloseCountGenome <- cbind(celluloseCountGenome, NCBI.Taxon.ID = NA)






# Get Other meta data info and save them in same dataset

fileName <- "P://Users personal data/phd materials/datas/imgdb/taxontable14673_20-nov-2014.xls"
genomesMetaInfo <- read.csv(fileName, header=TRUE, sep='\t')

i <- 1
#tempOrg <- data.frame(matrix(ncol=1, nrow=length(genomes)))

for ( org in genomes){
  #print(org)
  # find the genome in 'genomesMetaInfo'
  tempOrg <- genomesMetaInfo[genomesMetaInfo$Genome.Name == org,]
  #tempOrg[i,1] <- as.matrix(a[1]) 
  celluloseCountGenome$NCBI.Taxon.ID[i] <- tempOrg$NCBI.Taxon.ID[1]
  
  
  i <- i + 1
  
  
}

save(celluloseCountGenome, file = 'celluloseCountGenome2.Rda')
save(genomesMetaInfo, file = 'genomesMetaInfo2.Rda')


rm(list = ls())

pathname <- "P://Users personal data/phd materials/datas/imgdb/bacterialcellulose/"
setwd(pathname)
load('celluloseCountGenome2.Rda')
load('genomesMetaInfo2.Rda')

# combine the number of genes in pathways and cellulose 
fileName <- "P://Users personal data/phd materials/datas/kegg/geneCountInPathways.txt"
genesInPathways <- read.csv(fileName, header=TRUE, sep='\t')

# remove duplicate just sake of easy computation
genesInPathways <- genesInPathways[!duplicated(genesInPathways$organism),]

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

sink('ListOfPathwaysCoef2.txt')
coefOfPathways
sink()

X <- coef(summary(glmfit))
coefOfPathways <- cbind(coefOfPathways, V4 = NA)
coefOfPathways$V4[1:nrow(X)] <- X[,4]


