### R code from vignette source 'MAplots.Rnw'

###################################################
### code chunk number 1: loadPackage
###################################################
library(affyPLM)
options(width=60)


###################################################
### code chunk number 2: makeExpressionSet
###################################################
require(affydata)
data(Dilution)
eset.Dilution <- rma(Dilution)


###################################################
### code chunk number 3: FirstUsageShow (eval = FALSE)
###################################################
## par(mfrow=c(2,2))
## MAplot(eset.Dilution)


###################################################
### code chunk number 4: FirstUsageDo
###################################################
png("MAplotFirstUse.png",height=8,width=8,pointsize=10,res=300,units="in")
par(mfrow=c(2,2))
MAplot(eset.Dilution)
dev.off()


###################################################
### code chunk number 5: smoothScatterShow (eval = FALSE)
###################################################
## par(mfrow=c(2,2))
## MAplot(eset.Dilution,plot.method="smoothScatter")


###################################################
### code chunk number 6: smoothScatterDo
###################################################
##bitmap("MAplotsmoothScatter.png",height=8,width=8,pointsize=10,res=300)
png("MAplotsmoothScatter.png",height=8,width=8,pointsize=10,res=300,units="in")
par(mfrow=c(2,2))
MAplot(eset.Dilution,plot.method="smoothScatter",nrpoints=256)
dev.off()


###################################################
### code chunk number 7: RefUsageShow (eval = FALSE)
###################################################
## par(mfrow=c(2,2))
## MAplot(eset.Dilution,plot.method="smoothScatter",ref=1)


###################################################
### code chunk number 8: RefUsageDo
###################################################
png("MAplotsmoothScatterRef.png",height=8,width=8,pointsize=10,res=300,units="in")
par(mfrow=c(2,2))
MAplot(eset.Dilution,ref=1,plot.method="smoothScatter",nrpoints=256)
dev.off()


###################################################
### code chunk number 9: RefUsageWhichShow (eval = FALSE)
###################################################
## par(mfrow=c(2,1))
## MAplot(eset.Dilution,which=c(2,4),ref=1,plot.method="smoothScatter")


###################################################
### code chunk number 10: RefUsageWhichShow
###################################################
png("MAplotsmoothScatterWhich.png",height=8,width=8,pointsize=10,res=300,units="in")
par(mfrow=c(2,1))
MAplot(eset.Dilution,which=c(2,4),ref=1,plot.method="smoothScatter",nrpoints=256)
dev.off()


###################################################
### code chunk number 11: SubsetUsageWhichShow (eval = FALSE)
###################################################
## par(mfrow=c(2,1))
## MAplot(eset.Dilution,which=c(1,2),ref=c(1,2),plot.method="smoothScatter")


###################################################
### code chunk number 12: SubsetUsageWhichDo
###################################################
png("MAplotsmoothScatterSubset.png",height=8,width=8,pointsize=10,res=300,units="in")
par(mfrow=c(2,1))
MAplot(eset.Dilution,which=c(1,2),ref=c(1,2),plot.method="smoothScatter",nrpoints=256)
dev.off()


###################################################
### code chunk number 13: SubsetUsageWithNames
###################################################
MAplot(eset.Dilution,which=c("20A","20B"),ref=c("20A","20B"),plot.method="smoothScatter",nrpoints=256)


###################################################
### code chunk number 14: calcPA
###################################################
PA.calls <- mas5calls(Dilution)
Is.Present <- exprs(PA.calls) == "P"
Number.Present <- apply(Is.Present,1,sum)


###################################################
### code chunk number 15: addShow (eval = FALSE)
###################################################
## MAplot(eset.Dilution[Number.Present ==4,],show.statistics=FALSE,which=1,pch=20,cex=0.4,ylim=c(-0.8,0.8),xlim=c(2,15),add.loess=FALSE)
## MAplot(eset.Dilution[Number.Present ==0,],plot.method="add",col="red",which=1,pch=20,cex=0.4,show.statistics=FALSE,add.loess=FALSE)
## MAplot(eset.Dilution[Number.Present ==3,],plot.method="add",show.statistics=FALSE,which=1,col="green",pch=20,cex=0.4,add.loess=FALSE)
## MAplot(eset.Dilution[Number.Present ==2,],plot.method="add",show.statistics=FALSE,which=1,col="blue",pch=20,cex=0.4,add.loess=FALSE)
## MAplot(eset.Dilution[Number.Present ==1,],plot.method="add",show.statistics=FALSE,which=1,col="orange",pch=20,cex=0.4,add.loess=FALSE)


###################################################
### code chunk number 16: addDo
###################################################
png("MAplotadd.png",height=8,width=8,pointsize=10,res=300,units="in")
par(mfrow=c(1,1))
MAplot(eset.Dilution[Number.Present ==4,],show.statistics=FALSE,which=1,pch=20,cex=0.4,ylim=c(-0.8,0.8),xlim=c(2,15),add.loess=FALSE)
MAplot(eset.Dilution[Number.Present ==0,],plot.method="add",col="red",which=1,pch=20,cex=0.4,show.statistics=FALSE,add.loess=FALSE)
MAplot(eset.Dilution[Number.Present ==3,],plot.method="add",show.statistics=FALSE,which=1,col="green",pch=20,cex=0.4,add.loess=FALSE)
MAplot(eset.Dilution[Number.Present ==2,],plot.method="add",show.statistics=FALSE,which=1,col="blue",pch=20,cex=0.4,add.loess=FALSE)
MAplot(eset.Dilution[Number.Present ==1,],plot.method="add",show.statistics=FALSE,which=1,col="orange",pch=20,cex=0.4,add.loess=FALSE)
dev.off()


###################################################
### code chunk number 17: groupsShow (eval = FALSE)
###################################################
## MAplot(eset.Dilution,groups=c("Liver 20","Liver 20","Liver 10","Liver 10"),ref="Liver 10")


###################################################
### code chunk number 18: groupsDo
###################################################
png("MAplotgroups.png",height=8,width=8,pointsize=10,res=300,units="in")
MAplot(eset.Dilution,groups=c("Liver 20","Liver 20","Liver 10","Liver 10"),ref="Liver 10")
dev.off()


###################################################
### code chunk number 19: MAplots.Rnw:257-259
###################################################
## give ghostscript on Windows a few seconds to catch up
Sys.sleep(10)


