### R code from vignette source 'twoChannels.Rnw'

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
experimentName <- "DualChannelScreen"
dataPath <- system.file(experimentName, package="cellHTS2") 


###################################################
### code chunk number 4: readPlateList
###################################################
x <- readPlateList("Platelist.txt", name=experimentName, path=dataPath)


###################################################
### code chunk number 5: showX
###################################################
x


###################################################
### code chunk number 6: plateFileTable
###################################################
cellHTS2:::tableOutput(file.path(dataPath, "Platelist.txt"), "plate list", 
                       preName="twoChannels")


###################################################
### code chunk number 7: configure the data
###################################################
x <- configure(x, "Description.txt", "Plateconf.txt", "Screenlog.txt",
               path=dataPath) 


###################################################
### code chunk number 8: plateConfscreenLogTable
###################################################
cellHTS2:::tableOutputWithHeaderRows(file.path(dataPath, "Plateconf.txt"), 
                                     "plate configuration", selRows=NULL, 
                                     preName="twoChannels")

cellHTS2:::tableOutput(file.path(dataPath, "Screenlog.txt"), 
                       "screen log", selRows=1:2, preName="twoChannels")


###################################################
### code chunk number 9: twoChannels.Rnw:198-199
###################################################
table(wellAnno(x))


###################################################
### code chunk number 10: define controls
###################################################
## Define the controls for the different channels:
negControls <- vector("character", length=dim(Data(x))[3])

# channel 1 - gene A
negControls[1] <- "(?i)^geneA$" 
## case-insensitive and match the empty string at the beginning and 
## end of a line (to distinguish between "geneA" and "geneAB", for example.
## Although it is not a problem for the present well annotation)
 
# channel 2 - gene A and geneB
negControls[2] <- "(?i)^geneA$|^geneB$" 

posControls <- vector("character", length=dim(Data(x))[3])
# channel 1 - no controls
# channel 2 - geneC and geneD
posControls[2] <- "(?i)^geneC$|^geneD$"


###################################################
### code chunk number 11: writeReport1Show (eval = FALSE)
###################################################
## out <- writeReport(raw=x, outdir="raw",
##                    posControls=posControls, negControls=negControls)


###################################################
### code chunk number 12: writeReport1Do
###################################################
out <- writeReport(raw=x, force=TRUE, outdir="raw", 
   posControls=posControls, negControls=negControls)


###################################################
### code chunk number 13: browseReport1 (eval = FALSE)
###################################################
## if (interactive()) browseURL(out)


###################################################
### code chunk number 14: plateMedianChannels
###################################################
xn <- normalizePlates(x, scale="multiplicative", method="median", 
                      varianceAdjust="none")


###################################################
### code chunk number 15: set cut-off for R1
###################################################
ctoff <- quantile(Data(xn)[,,1], probs=0.05, na.rm=TRUE)


###################################################
### code chunk number 16: FvsRcorrected
###################################################
R <- Data(xn)[,,1]
F <- Data(xn)[,,2]

# Use the controls of R2 channel:
posC <- which(regexpr(posControls[2], as.character(wellAnno(x)), perl=TRUE)>0)
negC <- which(regexpr(negControls[2], as.character(wellAnno(x)), perl=TRUE)>0)

ylim <- range(F, na.rm=TRUE)
xlim <- range(R, na.rm=TRUE)

par(mfrow=c(1,ncol(F)), mai=c(1.15,1.15, 0.3,0.3))
for (r in 1:ncol(F)) {
#ind <- apply(cbind(R[,r],F[,r]), 1, function(z) any(is.na(z)))
#plot(R[!ind,r],F[!ind,r], col= densCols(cbind(R[!ind,r], F[!ind,r])), pch=16,cex=0.7,
#     xlab="R1 (log scale)", ylab="R2 (log scale)", log="xy", 
#    ylim=ylim, xlim=xlim)

plot(R[,r],F[,r], col= densCols(cbind(R[,r], F[,r])), pch=16,cex=0.7,
     xlab="R1 (log scale)", ylab="R2 (log scale)", log="xy", 
     ylim=ylim, xlim=xlim)

abline(v=ctoff, col="grey", lty=2, lwd=2)
points(R[posC,r],F[posC,r], col="red", pch=20, cex=0.9)
points(R[negC,r],F[negC,r], col="green", pch=20, cex=0.9)
ind <- which(Data(xn)[,r,1] <= ctoff)
points(R[ind,r],F[ind,r], col="grey", pch=20, cex=0.9)
#legend("topleft", col=c("grey", "red", "green"), legend=c("masked", 
#"positive controls","negative controls"),  bty="n", pch=20, cex=0.9)
 }


###################################################
### code chunk number 17: summarizeChannels
###################################################
xn1 <- summarizeChannels(xn, fun = function(r1, r2, 
                             thresh=quantile(r1, probs=0.05, na.rm=TRUE)) 
                         ifelse(r1>thresh, r2/r1, as.numeric(NA)))


###################################################
### code chunk number 18: show summarized object
###################################################
dim(Data(xn1))


###################################################
### code chunk number 19: redo plate correction
###################################################
xn1 <- normalizePlates(xn1, scale="multiplicative", log=TRUE, method="median", 
                       varianceAdjust="none") 


###################################################
### code chunk number 20: score and summarize replicates
###################################################
xsc <- scoreReplicates(xn1, sign="-", method="zscore") 
xsc <- summarizeReplicates(xsc, summary="mean") 


###################################################
### code chunk number 21: scores
###################################################
par(mfrow=c(1,2))
ylim <- quantile(Data(xsc), c(0.001, 0.999), na.rm=TRUE)
boxplot(Data(xsc) ~ wellAnno(xsc), col="lightblue", outline=FALSE, ylim=ylim)
imageScreen(xsc, zrange=c(-2,4))


###################################################
### code chunk number 22: RedefineControls
###################################################
## Define the controls for the normalized intensities (only one channel):
# For the single channel, the negative controls are geneA and geneB 
negControls <- "(?i)^geneA$|^geneB$" 
posControls <- "(?i)^geneC$|^geneD$"


###################################################
### code chunk number 23: report2Show (eval = FALSE)
###################################################
## setSettings(list(platelist=list(intensities=list(include=TRUE)),
##                  screenSummary=list(scores=list(range=c(-4,4)))))
## out <- writeReport(raw=x, normalized=xn1, scored=xsc, 
##                    outdir="logRatio", 
##                    map=TRUE, 
##                    posControls=posControls, negControls=negControls)


###################################################
### code chunk number 24: report2Do
###################################################
setSettings(list(platelist=list(intensities=list(include=TRUE)),
                 screenSummary=list(scores=list(range=c(-4,4)))))
out <- writeReport(raw=x, normalized=xn1, scored=xsc, 
                   force=TRUE, outdir="logRatio", 
                   map=TRUE,
                   posControls=posControls, negControls=negControls)


###################################################
### code chunk number 25: browse2 (eval = FALSE)
###################################################
## if (interactive()) browseURL(out)


###################################################
### code chunk number 26: savex
###################################################
save(xsc, file=paste(experimentName, ".rda", sep=""))


###################################################
### code chunk number 27: sessionInfo
###################################################
toLatex(sessionInfo())


