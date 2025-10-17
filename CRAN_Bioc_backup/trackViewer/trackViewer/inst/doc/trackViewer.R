## ----echo=FALSE, results="hide", warning=FALSE--------------------------------
suppressPackageStartupMessages({
    library(trackViewer)
    library(rtracklayer)
    library(Gviz)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(org.Hs.eg.db)
    library(VariantAnnotation)
  library(httr)
})
knitr::opts_chunk$set(warning=FALSE, message=FALSE)

## ----plotComp,echo=TRUE,fig.keep='none'---------------------------------------
library(Gviz)
library(rtracklayer)
library(trackViewer)
extdata <- system.file("extdata", package="trackViewer",
                       mustWork=TRUE)
gr <- GRanges("chr11", IRanges(122929275, 122930122), strand="-")
fox2 <- importScore(file.path(extdata, "fox2.bed"), format="BED",
                    ranges=gr)
fox2$dat <- coverageGR(fox2$dat)

viewTracks(trackList(fox2), gr=gr, autoOptimizeStyle=TRUE, newpage=FALSE)

dt <- DataTrack(range=fox2$dat[strand(fox2$dat)=="-"] , 
                genome="hg19", type="hist", name="fox2", 
                window=-1, chromosome="chr11", 
                fill.histogram="black", col.histogram="NA",
                background.title="white",
                col.frame="white", col.axis="black",
                col="black", col.title="black")
plotTracks(dt, from=122929275, to=122930122, strand="-")

## ----Gviz,echo=FALSE,fig.cap='Plot data with **Gviz** and **trackViewer**. Please note that **trackViewer** can generate similar figure as **Gviz** with several lines of simple codes.',fig.width=8,fig.height=3----
viewerStyle <- trackViewerStyle()
setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .13, .02, .02))
empty <- DataTrack(showAxis=FALSE, showTitle=FALSE, background.title="white")
plotTracks(list(empty, dt), from=122929275, to=122930122, strand="-")
pushViewport(viewport(0, .5, 1, .5, just=c(0, 0)))
viewTracks(trackList(fox2), viewerStyle=viewerStyle, 
                 gr=gr, autoOptimizeStyle=TRUE, newpage=FALSE)
popViewport()
grid.text(label="Gviz track", x=.3, y=.4)
grid.text(label="trackViewer track", x=.3, y=.9)

## ----lostcode,echo=TRUE,fig.keep='none'---------------------------------------
gr <- GRanges("chr1", IRanges(c(1, 6, 10), c(3, 6, 12)), score=c(3, 4, 1))
dt <- DataTrack(range=gr, data="score", type="hist")
plotTracks(dt, from=2, to=11)
tr <- new("track", dat=gr, type="data", format="BED")
viewTracks(trackList(tr), chromosome="chr1", start=2, end=11)

## ----GvizLost,echo=FALSE,fig.cap='Plot data with **Gviz** and **trackViewer**. Note that **trackViewer** is not only including more details but also showing all the data involved in the given range.',fig.width=8,fig.height=3----
plotTracks(list(empty, dt), from=2, to=11)
pushViewport(viewport(0, .5, 1, .5, just=c(0, 0)))
viewTracks(trackList(tr), viewerStyle=viewerStyle, 
           chromosome="chr1", start=2, end=11, 
           autoOptimizeStyle=TRUE, newpage=FALSE)
popViewport()
grid.text(label="Gviz track", x=.3, y=.4)
grid.text(label="trackViewer track", x=.3, y=.9)

## ----importData---------------------------------------------------------------
library(trackViewer)
extdata <- system.file("extdata", package="trackViewer",
                       mustWork=TRUE)
repA <- importScore(file.path(extdata, "cpsf160.repA_-.wig"),
                    file.path(extdata, "cpsf160.repA_+.wig"),
                    format="WIG")
## Because the wig file does not contain any strand info, 
## we need to set it manually.
strand(repA$dat) <- "-"
strand(repA$dat2) <- "+"

## ----coverage-----------------------------------------------------------------
fox2 <- importScore(file.path(extdata, "fox2.bed"), format="BED",
                    ranges=GRanges("chr11", IRanges(122830799, 123116707)))
dat <- coverageGR(fox2$dat)
## We can split the data by strand into two different track channels
## Here, we set the dat2 slot to save the negative strand info. 
 
fox2$dat <- dat[strand(dat)=="+"]
fox2$dat2 <- dat[strand(dat)=="-"]

## ----geneModel----------------------------------------------------------------
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
gr <- GRanges("chr11", IRanges(122929275, 122930122), strand="-")
trs <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg19.knownGene,
                         org.Hs.eg.db,
                         gr=gr)

## ----geneTrack----------------------------------------------------------------
entrezIDforFMR1 <- get("FMR1", org.Hs.egSYMBOL2EG)
theTrack <- geneTrack(entrezIDforFMR1,TxDb.Hsapiens.UCSC.hg19.knownGene)[[1]]

## ----viewTracks,fig.cap='plot data and annotation information along genomic coordinates',fig.width=8,fig.height=3----
viewerStyle <- trackViewerStyle()
setTrackViewerStyleParam(viewerStyle, "margin", c(.1, .05, .02, .02))
trackList <- trackList(repA, fox2, trs)
vp <- viewTracks(trackList, 
                 gr=gr, viewerStyle=viewerStyle, 
                 autoOptimizeStyle=TRUE)
addGuideLine(c(122929767, 122929969), vp=vp)
addArrowMark(list(x=122929650, 
                  y=2), # 2 means track 2 from the bottom.
             label="label",
             col="blue",
             vp=vp)

## ----browseTrack,fig.cap='interactive tracks',fig.width=6,fig.height=4--------
browseTracks(trackList, gr=gr)

## ----plotgeneTrack,fig.cap='Plot multiple genes in one track',fig.width=6,fig.height=4----
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
grW <- parse2GRanges("chr11:122,830,799-123,116,707")
ids <- getGeneIDsFromTxDb(grW, TxDb.Hsapiens.UCSC.hg19.knownGene)
symbols <- mget(ids, org.Hs.egSYMBOL)
genes <- geneTrack(ids, TxDb.Hsapiens.UCSC.hg19.knownGene, 
                   symbols, asList=FALSE)
optSty <- optimizeStyle(trackList(repA, fox2, genes), theme="safe")
trackListW <- optSty$tracks
viewerStyleW <- optSty$style
viewTracks(trackListW, gr=grW, viewerStyle=viewerStyleW)

## ----viewTracksOperator1,fig.cap='show data with operator "+"',fig.width=8,fig.height=4----
newtrack <- repA
## Must keep the same format for dat and dat2
newtrack <- parseWIG(newtrack, "chr11", 122929275, 122930122)
newtrack$dat2 <- newtrack$dat
newtrack$dat <- fox2$dat2
setTrackStyleParam(newtrack, "color", c("#4A95E0", "#CF5D6D"))
viewTracks(trackList(newtrack, trs), 
           gr=gr, viewerStyle=viewerStyle, operator="+")

## ----viewTracksOperator2,fig.cap='show data with operator "-"',fig.width=8,fig.height=4----
viewTracks(trackList(newtrack, trs), gr=gr, viewerStyle=viewerStyle, operator="-")

## ----viewTracksOperator3,fig.cap='show data with multiple operators', fig.width=8,fig.height=5----
viewTracks(trackList(newtrack, newtrack, newtrack, trs), gr=gr, viewerStyle=viewerStyle, operator=c(NA, "-", "+"))

## ----viewTracksOperator4,fig.cap='show data with operator "-"',fig.width=8,fig.height=4----
newtrack$dat <- GRoperator(newtrack$dat, newtrack$dat2, col="score", operator="-")
newtrack$dat2 <- GRanges()
viewTracks(trackList(newtrack, trs), gr=gr, viewerStyle=viewerStyle)

## -----------------------------------------------------------------------------
library(trackViewer)
features <- GRanges("chr1", IRanges(c(1, 501, 1001), 
                                    width=c(120, 400, 405),
                                    names=paste0("block", 1:3)),
                    fill = c("#FF8833", "#51C6E6", "#DFA32D"),
                    height = c(0.02, 0.05, 0.08))
SNP <- c(10, 100, 105, 108, 400, 410, 420, 600, 700, 805, 840, 1400, 1402)
sample.gr <- GRanges("chr1", IRanges(SNP, width=1, names=paste0("snp", SNP)),
                     color = sample.int(6, length(SNP), replace=TRUE),
                     score = sample.int(5, length(SNP), replace = TRUE))
lolliplot(sample.gr, features)


## ----fig.width=4.5,fig.height=3-----------------------------------------------
library(trackViewer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
methy <- import(system.file("extdata", "methy.bed", package="trackViewer"), "BED")
gr <- GRanges("chr22", IRanges(50968014, 50970514, names="TYMP"))
trs <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg19.knownGene,
                         org.Hs.eg.db,
                         gr=gr)
features <- c(range(trs[[1]]$dat), range(trs[[5]]$dat))
names(features) <- c(trs[[1]]$name, trs[[5]]$name)
features$fill <- c("lightblue", "mistyrose")
features$height <- c(.02, .04)
dandelion.plot(methy, features, ranges=gr, type="pin")

## -----------------------------------------------------------------------------
library(InteractionSet)
gi <- readRDS(system.file("extdata", "nij.chr6.51120000.53200000.gi.rds", package="trackViewer"))
head(gi)
range <- GRanges("chr6", IRanges(51120000, 53200000))
tr <- gi2track(gi)
ctcf <- readRDS(system.file("extdata", "ctcf.sample.rds",
                            package="trackViewer"))
## change the color
setTrackStyleParam(tr, "breaks", 
                   c(seq(from=0, to=50, by=10), 200))
setTrackStyleParam(tr, "color",
                   c("lightblue", "yellow", "red"))
viewTracks(trackList(ctcf, tr, heightDist=c(1, 3)), gr=range,
           autoOptimizeStyle = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  ideo <- loadIdeogram("hg38", chrom=c("chr1", "chr3", "chr4", 'chr21', "chr22"))

## ----echo=FALSE, results="hide"-----------------------------------------------
path <- system.file("extdata", "ideo.hg38.rds", package = "trackViewer")
ideo <- readRDS(path)

## ----eval=FALSE---------------------------------------------------------------
#  dataList <- trim(ideo)
#  dataList$score <- as.numeric(as.factor(dataList$gieStain))
#  dataList <- dataList[dataList$gieStain!="gneg"]
#  dataList <- GRangesList(dataList)
#  ideogramPlot(ideo, dataList,
#              layout=list("chr1", c("chr3", "chr22"),
#                          c("chr4", "chr21")))

## ----sessionInfo, results='asis'----------------------------------------------
sessionInfo()

