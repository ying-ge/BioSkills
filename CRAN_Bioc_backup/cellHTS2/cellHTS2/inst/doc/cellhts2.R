### R code from vignette source 'cellhts2.Rnw'

###################################################
### code chunk number 1: Ropts
###################################################
options(width=70)


###################################################
### code chunk number 2: install (eval = FALSE)
###################################################
## if (!requireNamespace("BiocManager", quietly=TRUE))
##     install.packages("BiocManager")
## BiocManager::install("cellHTS2")


###################################################
### code chunk number 3: setup1
###################################################
library("cellHTS2")


###################################################
### code chunk number 4: dataPath
###################################################
experimentName <- "KcViab"
dataPath <- system.file(experimentName, package="cellHTS2")


###################################################
### code chunk number 5: dirDataPath
###################################################
dataPath
rev(dir(dataPath))[1:12]


###################################################
### code chunk number 6: readPlateList
###################################################
x <- readPlateList("Platelist.txt",
                   name=experimentName,
                   path=dataPath)


###################################################
### code chunk number 7: showX
###################################################
x


###################################################
### code chunk number 8: plateFileTable
###################################################
cellHTS2:::tableOutput(file.path(dataPath, "Platelist.txt"), "plate list")
cellHTS2:::tableOutput(file.path(dataPath, names(intensityFiles(x))[1]),
              "signal intensity", header=FALSE)


###################################################
### code chunk number 9: see object state
###################################################
state(x)


###################################################
### code chunk number 10: writeReport
###################################################
out <- writeReport(raw=x, force = TRUE, outdir = "report-raw")


###################################################
### code chunk number 11: printout
###################################################
out


###################################################
### code chunk number 12: browseReport1 (eval = FALSE)
###################################################
## if (interactive()) browseURL(out)


###################################################
### code chunk number 13: annotatePlateRes
###################################################
x <- configure(x,
               descripFile="Description.txt",
               confFile="Plateconf.txt",
               logFile="Screenlog.txt",
               path=dataPath)


###################################################
### code chunk number 14: plateConfscreenLogTable
###################################################
cellHTS2:::tableOutputWithHeaderRows(file.path(dataPath, "Plateconf.txt"),
  "plate configuration", selRows=NULL)
cellHTS2:::tableOutput(file.path(dataPath, "Screenlog.txt"),
  "screen log", selRows=1:3)


###################################################
### code chunk number 15: cellhts2.Rnw:581-582
###################################################
table(wellAnno(x))


###################################################
### code chunk number 16: configurationplot
###################################################
configurationAsScreenPlot(x, legend=TRUE)


###################################################
### code chunk number 17: normalizePlateMedian
###################################################
xn <- normalizePlates(x,
                      scale="multiplicative",
                      log=FALSE,
                      method="median",
                      varianceAdjust="none")


###################################################
### code chunk number 18: compare cellHTs objects
###################################################
compare2cellHTS(x, xn)


###################################################
### code chunk number 19: score replicates
###################################################
xsc <- scoreReplicates(xn, sign="-", method="zscore")


###################################################
### code chunk number 20: summarize replicates
###################################################
xsc <- summarizeReplicates(xsc, summary="mean")


###################################################
### code chunk number 21: boxplotzscore
###################################################
scores <- Data(xsc)
ylim <- quantile(scores, c(0.001, 0.999), na.rm=TRUE)
boxplot(scores ~ wellAnno(x),
        col="lightblue", outline=FALSE, ylim=ylim)


###################################################
### code chunk number 22: callvalues
###################################################
y <- scores2calls(xsc, z0=1.5, lambda=2)
png("cellhts2-callvalues.png")
plot(Data(xsc), Data(y), col="blue", pch=".",
     xlab="z-scores", ylab="calls",
     main=expression(1/(1+e^{-lambda *(z-z[0])})))
dev.off()


###################################################
### code chunk number 23: callvaluesShow (eval = FALSE)
###################################################
## y <- scores2calls(xsc, z0=1.5, lambda=2)
## plot(Data(xsc), Data(y), col="blue", pch=".",
##      xlab="z-scores", ylab="calls",
##      main=expression(1/(1+e^{-lambda *(z-z[0])})))


###################################################
### code chunk number 24: geneIDs
###################################################
xsc <- annotate(xsc, geneIDFile="GeneIDs_Dm_HFA_1.1.txt",
                path=dataPath)


###################################################
### code chunk number 25: geneIDsTable
###################################################
cellHTS2:::tableOutput(file.path(dataPath, "GeneIDs_Dm_HFA_1.1.txt"),
     "gene ID", selRows = 3:6)


###################################################
### code chunk number 26: printxagain
###################################################
xsc


###################################################
### code chunk number 27: savex
###################################################
save(xsc, file=paste(experimentName, ".rda", sep=""))


###################################################
### code chunk number 28: writeReport2
###################################################
setSettings(list(plateList=list(reproducibility=list(include=TRUE, map=TRUE),
                                intensities=list(include=TRUE, map=TRUE)),
                 screenSummary=list(scores=list(range=c(-4, 8), map=TRUE))))
out = writeReport(raw=x, normalized=xn, scored=xsc,
                  force = TRUE, outdir = "report-normalized")


###################################################
### code chunk number 29: browseReport2 (eval = FALSE)
###################################################
## if (interactive()) browseURL(out)


###################################################
### code chunk number 30: imageScreen (eval = FALSE)
###################################################
## imageScreen(xsc, ar=1, zrange=c(-3,4))


###################################################
### code chunk number 31: exportData (eval = FALSE)
###################################################
## writeTab(xsc, file="Scores.txt")


###################################################
### code chunk number 32: exportOtherData (eval = FALSE)
###################################################
## # determine the ratio between each well and the plate median
## y <- array(as.numeric(NA), dim=dim(Data(x)))
## nrWell <- prod(pdim(x))
## nrPlate <- max(plate(x))
## for(p in 1:nrPlate)
## {
##     j <- (1:nrWell)+nrWell*(p-1)
##     samples <- wellAnno(x)[j]=="sample"
##     y[j, , ] <- apply(Data(x)[j, , , drop=FALSE], 2:3,
##                       function(w) w/median(w[samples],
##                                            na.rm=TRUE))
## }
## 
## y+signif(y, 3)
## out <- y[,,1]
## out <- cbind(fData(xsc), out)
## names(out) <- c(names(fData(xsc)),
## sprintf("Well/Median_r%d_ch%d", rep(1:dim(y)[2], dim(y)[3]),
## rep(1:dim(y)[3], each=dim(y)[2])))
## write.tabdel(out, file="WellMedianRatio.txt")


###################################################
### code chunk number 33: example for description file
###################################################
out <- templateDescriptionFile("template-Description.txt",
                               force=TRUE)
out
readLines(out)


###################################################
### code chunk number 34: old plateConfscreenLogTable
###################################################
cellHTS2:::tableOutput(file.path(dataPath, "old-Plateconf.txt"),
  "cellHTS package-specific plate configuration", selRows=1:28)
cellHTS2:::tableOutput(file.path(dataPath, "old-Screenlog.txt"),
  "cellHTS package-specific screen log", selRows=1:3)


###################################################
### code chunk number 35: Z score method (eval = FALSE)
###################################################
##  xZ <- normalizePlates(x, scale="additive", log=FALSE,
##          method="median", varianceAdjust="byPlate")


###################################################
### code chunk number 36: transfplots
###################################################
library(vsn)
myPlots=function(z, main, plotCol)
{
  z <- as.data.frame(z)
  colnames(z) <- paste0("Sample", seq_len(ncol(z)))
  gh <- ggplot2::ggplot(z, ggplot2::aes(x=Sample1))+
    ggplot2::geom_histogram(fill="darkblue")+
    ggplot2::ggtitle(main)
  gm <- vsn::meanSdPlot(as.matrix(z), plot=FALSE)$gg+
    ggplot2::ylim(c(0, quantile(abs(z[,2]-z[,1]), 0.95, na.rm=TRUE)))+
    ggplot2::theme(legend.key.size=unit(0.02, "npc"), legend.position="top")
  gq <- ggplot2::qplot(sample=z$Sample1)
  print(gh, vp=viewport(layout.pos.row=1, layout.pos.col=plotCol))
  print(gm, vp=viewport(layout.pos.row=2, layout.pos.col=plotCol))
  print(gq, vp=viewport(layout.pos.row=3, layout.pos.col=plotCol))
  return(NULL)
}
png("cellhts2-transfplots.png", width=400, height=600)
grid.newpage()
pushViewport(viewport(layout=grid.layout(3, 2)))
dv <- Data(xn)[,,1]
myPlots(dv, main="untransformed", plotCol=1)
xlog <- normalizePlates(x, scale="multiplicative", log=TRUE,
                        method="median", varianceAdjust="byExperiment")
dvlog <- Data(xlog)[,,1]
myPlots(dvlog, main="log2", plotCol=2)
dev.off()


###################################################
### code chunk number 37: sessionInfo
###################################################
toLatex(sessionInfo())


