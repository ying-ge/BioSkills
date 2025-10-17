### R code from vignette source 'miRNAtap.Rnw'

###################################################
### code chunk number 1: miRNAtap.Rnw:73-76 (eval = FALSE)
###################################################
## if (!requireNamespace("BiocManager", quietly=TRUE))
##     install.packages("BiocManager")
## BiocManager::install()


###################################################
### code chunk number 2: miRNAtap.Rnw:82-86 (eval = FALSE)
###################################################
## if (!requireNamespace("BiocManager", quietly=TRUE))
##     install.packages("BiocManager")
## BiocManager::install("miRNAtap")
## BiocManager::install("miRNAtap.db")


###################################################
### code chunk number 3: miRNAtap.Rnw:101-105 (eval = FALSE)
###################################################
## if (!requireNamespace("BiocManager", quietly=TRUE))
##     install.packages("BiocManager")
## BiocManager::install("topGO")
## BiocManager::install("org.Hs.eg.db")


###################################################
### code chunk number 4: miRNAtap.Rnw:110-113
###################################################
library(miRNAtap)
library(topGO)
library(org.Hs.eg.db)


###################################################
### code chunk number 5: miRNAtap.Rnw:119-122
###################################################
mir = 'miR-10b'
predictions = getPredictedTargets(mir, species = 'hsa',
                                    method = 'geom', min_src = 2)


###################################################
### code chunk number 6: miRNAtap.Rnw:127-128
###################################################
head(predictions)


###################################################
### code chunk number 7: miRNAtap.Rnw:137-140
###################################################
predictions_min = getPredictedTargets(mir, species = 'hsa',
                                    method = 'min', min_src = 2)
head(predictions_min)


###################################################
### code chunk number 8: miRNAtap.Rnw:149-151
###################################################
predictions_rat = getPredictedTargets(mir, species = 'rno',
                                    method = 'geom', min_src = 2)


###################################################
### code chunk number 9: miRNAtap.Rnw:158-166
###################################################
rankedGenes = predictions[,'rank_product']
selection = function(x) TRUE 
# we do not want to impose a cut off, instead we are using rank information
allGO2genes = annFUN.org(whichOnto='BP', feasibleGenes = NULL,
                mapping="org.Hs.eg.db", ID = "entrez")
GOdata =  new('topGOdata', ontology = 'BP', allGenes = rankedGenes, 
            annot = annFUN.GO2genes, GO2genes = allGO2genes, 
            geneSel = selection, nodeSize=10)


###################################################
### code chunk number 10: miRNAtap.Rnw:172-174
###################################################
results.ks = runTest(GOdata, algorithm = "classic", statistic = "ks")
results.ks


###################################################
### code chunk number 11: miRNAtap.Rnw:180-182
###################################################
allRes = GenTable(GOdata, KS = results.ks, orderBy = "KS", topNodes = 20)
allRes[,c('GO.ID','Term','KS')]


###################################################
### code chunk number 12: miRNAtap.Rnw:193-194
###################################################
toLatex(sessionInfo())


