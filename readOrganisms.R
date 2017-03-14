# Author: Syeda Sakira Hassan
# Date: 05.08.2014
# ChangeLog Dates: 25.11.2014

# Description: 

rm(list = ls())
library(KEGGREST)


## read the organism file & list them
fileName <- "P://Users personal data/phd materials/datas/kegg/organism.txt"
organisms <- read.csv(fileName, header=FALSE, sep='\t')
head(organisms)
nrow(organisms) # The number of organisms is 3175 / 4054 ( new update)

## Find the organisms with Bacteria Family

# Find the index
bacteriaCol <- grep("Bacteria", organisms$V4)
length(bacteriaCol) # The number of bacterial organisms is 2731 / 3526 (new update)

# Find the organisms
bacterialOrganisms <- organisms[bacteriaCol,]
head(bacterialOrganisms)
tail(bacterialOrganisms)

## Get the list of genes for each bacterial organism

n <- nrow(bacterialOrganisms)
n

#allGeneList <- data.frame()
allGeneListFile <- "P://Users personal data/phd materials/datas/kegg/all_gene_list.txt"
file.create(allGeneListFile)

# KEGG organism code
orgCodeList <- bacterialOrganisms$V2

########################################
for( i in 1:n){
	orgCode <- orgCodeList[i]
	geneList <- keggList(orgCode)
	write.table(geneList, file = allGeneListFile, sep="\t", append=TRUE, col.names = FALSE)

	#allGeneList <- rbind(allGeneList,geneList)
}

###############################################


#### Save all gene list in a directory
allGeneListDir <- "P://Users personal data/phd materials/datas/kegg/all_gene_list/"
for( i in 1:n){
	orgCode <- orgCodeList[i]
	geneList <- keggList(orgCode)
	allGeneListFile <- paste(allGeneListDir,orgCode,'.txt', sep="")

	write.table(geneList, file = allGeneListFile, sep="\t", append=TRUE, col.names = FALSE)

	
}




# Save all pathways in pathways directory
allPathwayListDir <- "P://Users personal data/phd materials/datas/kegg/pathways/"
for( i in 1:n){
	orgCode <- orgCodeList[i]
	pathwayList <- tryCatch({
		keggList("pathway", orgCode)
		},
		error = function(err) {
 
  			# error handler picks up where error was generated
  			# print(paste("MY_ERROR:  ",err))
			print( paste(err, "\nThe pathways for the organism ", orgCode, " was not found in KEGG database."))
  		}
	) # End tryCatch
	allPathwayListFile <- paste(allPathwayListDir,orgCode,'.txt', sep="")

	write.table(pathwayList, file = allPathwayListFile, sep="\t", append=TRUE, col.names = FALSE)

	
}


# Save all genes link to pathways in genes_pathways directory
allGenePathwayListDir <- "P://Users personal data/phd materials/datas/kegg/genes_pathways/"
for( i in 1:n){
	orgCode <- orgCodeList[i]
	genePathwayList <- tryCatch({
		keggLink( "pathway", orgCode )
		},
		error = function(err) {
 
  			# error handler picks up where error was generated
  			# print(paste("MY_ERROR:  ",err))
			print( paste(err, "\nThe genes link to pathways for the organism ", orgCode, " was not found in KEGG database."))
  		}
	) # End tryCatch
	allGenePathwayListFile <- paste(allGenePathwayListDir,orgCode,'.txt', sep="")

	write.table(cbind(genePathwayList), file = allGenePathwayListFile, sep="\t",  col.names = FALSE, row.names = TRUE)

	
}




