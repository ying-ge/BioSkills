### R code from vignette source 'Overview.Rnw'

###################################################
### code chunk number 1: Overview.Rnw:70-71
###################################################
library(pathifier) 


###################################################
### code chunk number 2: Overview.Rnw:77-78
###################################################
data(Sheffer)


###################################################
### code chunk number 3: Overview.Rnw:84-85
###################################################
data(KEGG)


###################################################
### code chunk number 4: Overview.Rnw:91-94
###################################################
PDS<-quantify_pathways_deregulation(sheffer$data, sheffer$allgenes,
  kegg$gs, kegg$pathwaynames, sheffer$normals, attempts = 100,
  min_exp=sheffer$minexp, min_std=sheffer$minstd)


###################################################
### code chunk number 5: Overview.Rnw:100-105
###################################################
x<-NULL
x$normals<-PDS$scores$MISMATCH_REPAIR[sheffer$normals]
x$tumors<-PDS$scores$MISMATCH_REPAIR[!sheffer$normals]
boxplot(x)
boxplot(x,ylab="score")


###################################################
### code chunk number 6: Overview.Rnw:109-110
###################################################
as.character(sheffer$samples[PDS$scores$REGULATION_OF_AUTOPHAGY>0.8])


