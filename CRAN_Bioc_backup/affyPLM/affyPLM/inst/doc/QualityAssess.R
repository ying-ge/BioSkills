### R code from vignette source 'QualityAssess.Rnw'

###################################################
### code chunk number 1: loadPackage
###################################################
library(affyPLM)
options(width=40)


###################################################
### code chunk number 2: loadData
###################################################
require(affydata)
data(Dilution)   # an example dataset provided by the affydata package
#FIXME:remove the next line
Dilution = updateObject(Dilution)
Pset <- fitPLM(Dilution)


###################################################
### code chunk number 3: weightsImageShow (eval = FALSE)
###################################################
## image(Pset,which=2)


###################################################
### code chunk number 4: weightsImageDo
###################################################
png("Quality-weightimage1a.png",height=4,width=4,pointsize=10,res=300,units="in")
par(mar=c(2.0,2.1,1.6,1.1),oma=c(1,1,0,0))
image(Pset,which=2)
dev.off()


###################################################
### code chunk number 5: weightscolorImageShow (eval = FALSE)
###################################################
## image(Pset,which=2,col=gray(0:25/25),add.legend=TRUE)
## image(Pset,which=2,col=gray(25:0/25),add.legend=TRUE)


###################################################
### code chunk number 6: weightscolorImageDo
###################################################
png("Quality-weightimage2a.png",height=8,width=8,res=300,units="in")
par(mar=c(2.0,2.1,1.6,1.1),oma=c(1,1,0,0))
image(Pset,which=2,col=gray(0:25/25),add.legend=T)
dev.off()

png("Quality-weightimage2b.png",height=8,width=8,res=300,units="in")
par(mar=c(2.0,2.1,1.6,1.1),oma=c(1,1,0,0))
image(Pset,which=2,col=gray(25:0/25),add.legend=T)
dev.off()


###################################################
### code chunk number 7: residualImageShow (eval = FALSE)
###################################################
## image(Pset,which=2, type="resids")
## image(Pset,which=2, type="pos.resids")
## image(Pset,which=2, type="neg.resids")
## image(Pset,which=2, type="sign.resids")


###################################################
### code chunk number 8: residualImageDo
###################################################
png("Quality-residualimages1.png",height=4,width=4,res=300,units="in")
par(mar=c(2.0,2.1,1.6,1.1),oma=c(1,1,0,0))
image(Pset,which=2, type="resids")
dev.off()
png("Quality-residualimages2.png",height=4,width=4,res=300,units="in")
par(mar=c(2.0,2.1,1.6,1.1),oma=c(1,1,0,0))
image(Pset,which=2, type="pos.resids")
dev.off()
png("Quality-residualimages3.png",height=4,width=4,res=300,units="in")
par(mar=c(2.0,2.1,1.6,1.1),oma=c(1,1,0,0))
image(Pset,which=2, type="neg.resids")
dev.off()
png("Quality-residualimages4.png",height=4,width=4,res=300,units="in")
par(mar=c(2.0,2.1,1.6,1.1),oma=c(1,1,0,0))
image(Pset,which=2, type="sign.resids")
dev.off()


###################################################
### code chunk number 9: residualcolorImageShow (eval = FALSE)
###################################################
## image(Pset,which=2,type="resids",col=pseudoPalette(low="darkgreen",high="magenta",mid="lightgrey"),add.legend=TRUE)
## image(Pset,which=2,type="pos.resids",col=pseudoPalette(low="yellow",high="darkblue"),add.legend=TRUE)


###################################################
### code chunk number 10: residualcolorImageDo
###################################################
png("Quality-residualimages5.png",height=8,width=8,res=300,units="in")
par(mar=c(2.0,2.1,1.6,1.1),oma=c(1,1,0,0))
image(Pset,which=2,type="resids",col=pseudoPalette(low="darkgreen",high="magenta",mid="lightgrey"),add.legend=TRUE)
dev.off()
png("Quality-residualimages6.png",height=8,width=8,res=300,units="in")
par(mar=c(2.0,2.1,1.6,1.1),oma=c(1,1,0,0))
image(Pset,which=2,type="pos.resids",col=pseudoPalette(low="yellow",high="darkblue"),add.legend=TRUE)
dev.off()


###################################################
### code chunk number 11: RLEShow (eval = FALSE)
###################################################
## RLE(Pset,main="RLE for Dilution dataset")


###################################################
### code chunk number 12: RLEDo
###################################################
png("Quality-RLE.png",height=4,width=4,pointsize=10,res=300,units="in")
par(mar=c(2.0,2.1,1.6,1.1),oma=c(1,1,0,0))
RLE(Pset,main="RLE for Dilution dataset")
dev.off()


###################################################
### code chunk number 13: rleStat
###################################################
RLE(Pset,type="stats")


###################################################
### code chunk number 14: NUSEShow (eval = FALSE)
###################################################
## NUSE(Pset,main="NUSE for Dilution dataset")


###################################################
### code chunk number 15: NUSEDo
###################################################
png("Quality-NUSE.png",height=4,width=4,pointsize=10,res=300,units="in")
par(mar=c(2.0,2.1,1.6,1.1),oma=c(1,1,0,0))
NUSE(Pset,main="NUSE for Dilution dataset")
dev.off()


###################################################
### code chunk number 16: nuseStat
###################################################
 NUSE(Pset,type="stats")


###################################################
### code chunk number 17: QualityAssess.Rnw:250-252
###################################################
## give ghostscript on Windows a few seconds to catch up (UGLY HACK)
Sys.sleep(10)


