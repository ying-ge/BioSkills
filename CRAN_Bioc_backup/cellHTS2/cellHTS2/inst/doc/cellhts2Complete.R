### R code from vignette source 'cellhts2Complete.Rnw'

###################################################
### code chunk number 1: Ropts
###################################################
options(width=70)


###################################################
### code chunk number 2: setup1
###################################################
library("cellHTS2")


###################################################
### code chunk number 3: setup2
###################################################
## working path:
workPath <- getwd()

## check if bib file exists
if (!("cellhts.bib" %in% dir()) )
  system(sprintf("cp %s/cellhts.bib .", system.file("doc", package="cellHTS2")))

## for debugging:
options(error=recover)
## for software development, when we do not want to install
## the package after each minor change:
##   for(f in dir("~/huber/projects/Rpacks/cellHTS2/R", full.names=TRUE, pattern=".R$"))source(f)


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
### code chunk number 10: writeReport1Show (eval = FALSE)
###################################################
## out <- writeReport(raw=x)


###################################################
### code chunk number 11: writeReport1Do
###################################################
out <- writeReport(raw=x, force=TRUE, outdir=tempdir())


###################################################
### code chunk number 12: printout
###################################################
out


###################################################
### code chunk number 13: browseReport1
###################################################
if (interactive()) browseURL(out)


###################################################
### code chunk number 14: annotatePlateRes
###################################################
x <- configure(x,
               descripFile="Description.txt",
               confFile="Plateconf.txt",
               logFile="Screenlog.txt",
               path=dataPath)


###################################################
### code chunk number 15: plateConfscreenLogTable
###################################################
cellHTS2:::tableOutputWithHeaderRows(file.path(dataPath, "Plateconf.txt"),
                                     "plate configuration", selRows=NULL)
cellHTS2:::tableOutput(file.path(dataPath, "Screenlog.txt"),
                       "screen log", selRows=1:3)


###################################################
### code chunk number 16: cellhts2Complete.Rnw:573-574
###################################################
table(wellAnno(x))


###################################################
### code chunk number 17: configurationplot
###################################################
png("cellhts2Complete-configurationplot.png", width=324, height=324)
configurationAsScreenPlot(x)
dev.off()


###################################################
### code chunk number 18: configurationplotShow (eval = FALSE)
###################################################
## configurationAsScreenPlot(x)


###################################################
### code chunk number 19: normalizePlateMedian
###################################################
xn <- normalizePlates(x,
                      scale="multiplicative",
                      log=FALSE,
                      method="median",
                      varianceAdjust="none")


###################################################
### code chunk number 20: compare cellHTs objects
###################################################
compare2cellHTS(x, xn)


###################################################
### code chunk number 21: score replicates
###################################################
xsc <- scoreReplicates(xn, sign="-", method="zscore")


###################################################
### code chunk number 22: summarize replicates
###################################################
xsc <- summarizeReplicates(xsc, summary="mean")


###################################################
### code chunk number 23: boxplotzscore
###################################################
scores <- Data(xsc)
ylim <- quantile(scores, c(0.001, 0.999), na.rm=TRUE)
boxplot(scores ~ wellAnno(x), col="lightblue", outline=FALSE,
        ylim=ylim)


###################################################
### code chunk number 24: callvalues
###################################################
y <- scores2calls(xsc, z0=1.5, lambda=2)
png("cellhts2Complete-callvalues.png")
plot(Data(xsc), Data(y), col="blue", pch=".",
     xlab="z-scores", ylab="calls",
     main=expression(1/(1+e^{-lambda *(z-z[0])})))
dev.off()


###################################################
### code chunk number 25: callvaluesShow (eval = FALSE)
###################################################
## y <- scores2calls(xsc, z0=1.5, lambda=2)
## plot(Data(xsc), Data(y), col="blue", pch=".",
##      xlab="z-scores", ylab="calls",
##      main=expression(1/(1+e^{-lambda *(z-z[0])})))


###################################################
### code chunk number 26: geneIDs
###################################################
xsc <- annotate(xsc, geneIDFile="GeneIDs_Dm_HFA_1.1.txt",
                path=dataPath)


###################################################
### code chunk number 27: geneIDsTable
###################################################
cellHTS2:::tableOutput(file.path(dataPath, "GeneIDs_Dm_HFA_1.1.txt"),
                       "gene ID", selRows = 3:6)


###################################################
### code chunk number 28: bdgpbiomart1
###################################################
data("bdgpbiomart")
fData(xsc) <- bdgpbiomart
fvarMetadata(xsc)[names(bdgpbiomart), "labelDescription"] <-
  sapply(names(bdgpbiomart),
    function(i) sub("_", " ", i)
)


###################################################
### code chunk number 29: get path for Rnw files (eval = FALSE)
###################################################
## rnwPath <- system.file("doc/Rnw", package="cellHTS2")
## setwd(rnwPath)
## system(sprintf("cp biomart.tex %s", workPath))
## setwd(workPath)


###################################################
### code chunk number 30: runBiomart (eval = FALSE)
###################################################
## cat("Weaving the biomart vignette part...")
## setwd(rnwPath)
## Sweave("biomart.Rnw")
## setwd(workPath)
## stop()


###################################################
### code chunk number 31: printxagain
###################################################
xsc


###################################################
### code chunk number 32: savex
###################################################
save(xsc, file=paste(experimentName, ".rda", sep=""))


###################################################
### code chunk number 33: writeReport2 (eval = FALSE)
###################################################
## setSettings(list(plateList=list(reproducibility=list(include=TRUE, map=TRUE),
##                                 intensities=list(include=TRUE, map=TRUE)),
##                  screenSummary=list(scores=list(range=c(-4, 8), map=TRUE))))
## out <- writeReport(raw=x, normalized=xn, scored=xsc,
##                    force=TRUE)


###################################################
### code chunk number 34: browseReport2 (eval = FALSE)
###################################################
## if (interactive()) browseURL(out)


###################################################
### code chunk number 35: imageScreen (eval = FALSE)
###################################################
## imageScreen(xsc, ar=1, zrange=c(-3,4))


###################################################
### code chunk number 36: exportData (eval = FALSE)
###################################################
## writeTab(xsc, file="Scores.txt")


###################################################
### code chunk number 37: exportOtherData (eval = FALSE)
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
## y <- signif(y, 3)
## out <- y[,,1]
## out <- cbind(fData(xsc), out)
## names(out) = c(names(fData(xsc)),
## sprintf("Well/Median_r%d_ch%d", rep(1:dim(y)[2], dim(y)[3]),
## rep(1:dim(y)[3], each=dim(y)[2])))
## write.tabdel(out, file="WellMedianRatio.txt")


###################################################
### code chunk number 38: category
###################################################
library("Category")


###################################################
### code chunk number 39: obsolete GO ids
###################################################
obsolete <- c("GO:0005489", "GO:0015997", "GO:0045034", "GO:0005660",
              "GO:0006118", "GO:0006512", "GO:0045045", "GO:0006125",
              "GO:0043072", "GO:0006100", "GO:0048740")


###################################################
### code chunk number 40: cat1
###################################################
scores <- as.vector(Data(xsc))
names(scores) <- geneAnno(xsc)
sel <- !is.na(scores) &
       (!is.na(fData(xsc)$go_biological_process_id))
goids <- strsplit(fData(xsc)$go_biological_process_id[sel],
                  ", ")
goids <- lapply(goids, function(x) x[!(x %in% obsolete)])
genes <- rep(geneAnno(xsc)[sel], listLen(goids))
categs <- cateGOry(genes, unlist(goids, use.names=FALSE))


###################################################
### code chunk number 41: cat2
###################################################
nrMem <- rowSums(categs) # number of genes per category
remGO <- which(nrMem < 3 | nrMem > 1000)
categs <- categs[-remGO,,drop=FALSE]
##  see if there are genes that don't belong to any category
## after applying the filter
nrMem <- rowSums(t(categs))
rem <- which(nrMem==0)
if(length(rem)!=0) categs <- categs[,-rem, drop=FALSE]


###################################################
### code chunk number 42: cat3
###################################################
stats <- scores[ sel & (names(scores) %in% colnames(categs)) ]


###################################################
### code chunk number 43: handle replicates
###################################################
## handle duplicated genes in stats:
isDup <- duplicated(names(stats))
table(isDup)
dupNames <- names(stats)[isDup]
sp <- stats[names(stats) %in% dupNames]
sp <- split(sp, names(sp))
table(sapply(sp, length))
aux <- stats[!isDup]
aux[names(sp)] <- sapply(sp, max)
stats <- aux
rm(aux)


###################################################
### code chunk number 44: arrange probes
###################################################
m <- match(colnames(categs), names(stats))
stats <- stats[m]
stopifnot(colnames(categs)==names(stats))


###################################################
### code chunk number 45: cat6
###################################################
acMean <- applyByCategory(stats, categs)
acTtest <- applyByCategory(stats, categs,
                           FUN=function(v)
                           t.test(v, stats)$p.value)
acNum <- applyByCategory(stats, categs, FUN=length)
isEnriched <- (acTtest<=1e-3) & (acMean>0.5)


###################################################
### code chunk number 46: volcano
###################################################
png("cellhts2Complete-volcano.png", width=180, height=180)
par(mai=c(0.9,0.9,0.1,0.1))
px <- cbind(acMean, -log10(acTtest))
plot(px, main='', xlab=expression(z[mean]),
     ylab=expression(-log[10]~p), pch=".", col="black")
points(px[isEnriched, ], pch=16, col="red", cex=0.7)
dev.off()
stopifnot(identical(names(acMean), names(acTtest)),
          identical(names(acMean), names(acNum)))


###################################################
### code chunk number 47: enrichedGoCateg
###################################################
enrichedGOCateg <- names(which(isEnriched))
require("GO.db")
res <- data.frame(
   "$n$" = acNum[isEnriched],
    "$z_{\\mbox{\\scriptsize mean}}$" = signif(acMean[isEnriched],2),
    "$p$" = signif(acTtest[isEnriched],2),
    "GOID" = I(enrichedGOCateg),
    "Ontology" = I(sapply(enrichedGOCateg, function(x) Ontology(get(x, GOTERM)))),
    "description" = I(sapply(enrichedGOCateg, function(x) Term(get(x, GOTERM)))),
    check.names=FALSE)

mt <- match(res$Ontology, c("CC", "BP", "MF"))
stopifnot(!any(is.na(mt)))
res <- res[order(mt, res$"$p$"), ]

cellHTS2:::dataframeOutput(res, header=TRUE,
  caption=sprintf("Top %d Gene Ontology categories with respect to $z$-score.", nrow(res)),
  label="enrichedGoCateg", gotable=TRUE)


###################################################
### code chunk number 48: load file with previous analysis
###################################################
data2003 <- read.table(file.path(dataPath, "Analysis2003.txt"),
                       header=TRUE, as.is=TRUE, sep="\t")


###################################################
### code chunk number 49: add the current scored values
###################################################
i <- data2003$Position + 384*(data2003$Plate-1)
data2003$ourScore <- as.vector(Data(xsc))[i]


###################################################
### code chunk number 50: scoresComparison
###################################################
png("cellhts2Complete-scoresComparison.png", width=324, height=324)
par(mfrow=c(7,9), mai=c(0,0,0,0))
for(i in 1:max(data2003$Plate))
{
   sel <- (data2003$Plate==i)
   plot(data2003$ourScore[sel], data2003$Score[sel], pch=19, cex=0.6)
}
dev.off()


###################################################
### code chunk number 51: example for description file
###################################################
out <- templateDescriptionFile("template-Description.txt",
                               force=TRUE)
out
readLines(out)


###################################################
### code chunk number 52: old plateConfscreenLogTable
###################################################
cellHTS2:::tableOutput(file.path(dataPath, "old-Plateconf.txt"),
                       "cellHTS package-specific plate configuration", selRows=1:28)
cellHTS2:::tableOutput(file.path(dataPath, "old-Screenlog.txt"),
                       "cellHTS package-specific screen log", selRows=1:3)


###################################################
### code chunk number 53: Z score method (eval = FALSE)
###################################################
## xZ <- normalizePlates(x, scale="additive", log=FALSE,
##                       method="median",
##                       varianceAdjust="byPlate")


###################################################
### code chunk number 54: transfplots
###################################################
library(vsn)
myPlots=function(z, main, plotCol)
{
  z <- as.data.frame(z)
  colnames(z) <- paste0("Sample", seq_len(ncol(z)))
  gh <- ggplot2::ggplot(z, ggplot2::aes(x=Sample1))+
    ggplot2::geom_histogram(fill="darkblue")+
    ggplot2::ggtitle(main)
  gm <- meanSdPlot(as.matrix(z), plot=FALSE)$gg+
    ggplot2::ylim(c(0, quantile(abs(z[,2]-z[,1]), 0.95, na.rm=TRUE)))+
    ggplot2::theme(legend.key.size=unit(0.02, "npc"), legend.position="top")
  gq <- ggplot2::qplot(sample=z$Sample1)
  print(gh, vp=viewport(layout.pos.row=1, layout.pos.col=plotCol))
  print(gm, vp=viewport(layout.pos.row=2, layout.pos.col=plotCol))
  print(gq, vp=viewport(layout.pos.row=3, layout.pos.col=plotCol))
}
png("cellhts2Complete-transfplots.png", width=400, height=600)
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
### code chunk number 55: sessionInfo
###################################################
toLatex(sessionInfo())


