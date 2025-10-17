### R code from vignette source 'twoWay.Rnw'

###################################################
### code chunk number 1: setup1
###################################################
library("cellHTS2")


###################################################
### code chunk number 2: setup2
###################################################
## for debugging:
options(error=recover)


###################################################
### code chunk number 3: dataPath
###################################################
experimentName <- "TwoWayAssay"
dataPath <- system.file(experimentName, package="cellHTS2") 


###################################################
### code chunk number 4: source import function
###################################################
source(file.path(dataPath, "importData.R"))


###################################################
### code chunk number 5: readPlateData
###################################################
x <- readPlateList("Platelist.txt", name=experimentName,
                   importFun=importData, path=dataPath)


###################################################
### code chunk number 6: showX
###################################################
x


###################################################
### code chunk number 7: plateFileTable
###################################################
cellHTS2:::tableOutput(file.path(dataPath, "Platelist.txt"), selRows=1, 
  "plate list", preName="twoWay")


###################################################
### code chunk number 8: configure the data
###################################################
x <- configure(x, descripFile = "Description.txt",
               confFile="Plateconf.txt", path=dataPath)


###################################################
### code chunk number 9: well annottaion
###################################################
table(wellAnno(x))


###################################################
### code chunk number 10: annotate the data
###################################################
x <- annotate(x, geneIDFile="GeneIDs.txt", path=dataPath)


###################################################
### code chunk number 11: plateConfscreenLogTable
###################################################
cellHTS2:::tableOutputWithHeaderRows(file.path(dataPath, "Plateconf.txt"), 
                                     "plate configuration", selRows=NULL, preName="twoWay")


###################################################
### code chunk number 12: geneIDsTable
###################################################
cellHTS2:::tableOutput(file.path(dataPath, "GeneIDs.txt"), "gene ID", 
                       selRows = 3:6, preName="twoWay")


###################################################
### code chunk number 13: define controls
###################################################
negCtr <- "(?i)^GFP$|^mock$"
posCtr <- list(act = "(?i)^AATK$|^ATTK$", inh = "(?i)^MAP2K6$")


###################################################
### code chunk number 14: writeReport1Show (eval = FALSE)
###################################################
## setSettings(list(platelist=list(intensities=list(range=c(300, 4000), 
##                                                  include=TRUE))))
## out <- writeReport(raw=x, outdir="2Wraw",
##                    posControls=posCtr, negControls=negCtr)


###################################################
### code chunk number 15: writeReport1Do
###################################################
setSettings(list(platelist=list(intensities=list(range=c(300, 4000), 
                                                 include=TRUE))))
out <- writeReport(raw=x, force=TRUE, outdir="2Wraw", 
                   posControls=posCtr, negControls=negCtr) 


###################################################
### code chunk number 16: browseReport1 (eval = FALSE)
###################################################
## browseURL(out)


###################################################
### code chunk number 17: normalization
###################################################
xn <- normalizePlates(x, scale="multiplicative", log=TRUE, method ="negatives", 
                      varianceAdjust="none", negControls = negCtr)


###################################################
### code chunk number 18: get the data as an array
###################################################
xnorm <- Data(xn)
dim(xnorm)


###################################################
### code chunk number 19: score and summarize replicates
###################################################
xsc <- scoreReplicates(xn, sign="+", method="zscore") 
xsc <- summarizeReplicates(xsc, summary="mean") 


###################################################
### code chunk number 20: boxplotzscore
###################################################
ylim <- quantile(Data(xsc), c(0.001, 0.999), na.rm=TRUE)
wa <- factor(as.character(wellAnno(xsc)), exclude="empty")   # to exclude "empty" wells
boxplot(Data(xsc) ~ wa, col="lightblue", main="scores", outline=FALSE, ylim=ylim, xaxt="n")
lab <- unique(plateConf(xsc)$Content)
lab <- lab[match(levels(wa), tolower(lab))]
axis(1, at=c(1:nlevels(wa)), labels=lab)


###################################################
### code chunk number 21: report2Show (eval = FALSE)
###################################################
## setSettings(list(platelist=list(intensities=list(range=c(-1, 1), include=TRUE)),
##                  screenSummary=list(scores=list(range=c(-2,3)))))
## out <- writeReport(raw=x, normalized=xn, scored=xsc, 
##                    outdir="2Wnormalized", posControls=posCtr, negControls=negCtr)


###################################################
### code chunk number 22: report2Do
###################################################
setSettings(list(platelist=list(intensities=list(range=c(-1, 1), include=TRUE)),
                 screenSummary=list(scores=list(range=c(-2,3)))))
out <- writeReport(raw=x, normalized=xn, scored=xsc, 
                   outdir="2Wnormalized", posControls=posCtr, negControls=negCtr, 
                   force=TRUE)


###################################################
### code chunk number 23: browse2 (eval = FALSE)
###################################################
## browseURL(out)


###################################################
### code chunk number 24: savex
###################################################
save(xsc, file=paste(experimentName, ".rda", sep=""))


###################################################
### code chunk number 25: sessionInfo
###################################################
toLatex(sessionInfo())


