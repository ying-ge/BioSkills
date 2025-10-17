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

## -----------------------------------------------------------------------------
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

fox2 <- importScore(file.path(extdata, "fox2.bed"), format="BED",
                    ranges=GRanges("chr11", IRanges(122830799, 123116707)))
dat <- coverageGR(fox2$dat)
## We can split the data by strand into two different track channels
## Here, we set the dat2 slot to save the negative strand info. 
 
fox2$dat <- dat[strand(dat)=="+"]
fox2$dat2 <- dat[strand(dat)=="-"]

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
gr <- GRanges("chr11", IRanges(122929275, 122930122), strand="-")
trs <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg19.knownGene,
                         org.Hs.eg.db,
                         gr=gr)


## ----optSty-------------------------------------------------------------------
optSty <- optimizeStyle(trackList(repA, fox2, trs))
trackList <- optSty$tracks
viewerStyle <- optSty$style

## ----browseTrack, fig.width=6,fig.height=4------------------------------------
browseTracks(trackList, gr=gr)

## ----viewTracksXaxis,fig.cap='plot data with x-scale',fig.width=8,fig.height=3----
setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .05, .01, .01))
setTrackXscaleParam(trackList[[1]], "draw", TRUE)
setTrackXscaleParam(trackList[[1]], "gp", list(cex=0.8))
viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
setTrackXscaleParam(trackList[[1]], attr="position", 
                    value=list(x=122929700, y=3, label=200))
viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)

## ----viewTracksYaxis,fig.cap='plot data with y-axis in right side',fig.width=8,fig.height=3----
setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .05, .01, .05))
for(i in 1:2){
    setTrackYaxisParam(trackList[[i]], "main", FALSE)
}
## Adjust the limit of y-axis
setTrackStyleParam(trackList[[1]], "ylim", c(0, 25))
setTrackStyleParam(trackList[[2]], "ylim", c(-25, 0))
viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)

## ----viewTracksYaxisCol,fig.cap='plot data with adjusted track color',fig.width=8,fig.height=3----
## change the 
setTrackYaxisParam(trackList[[1]], "gp", list(cex=.8, col="green"))
setTrackYaxisParam(trackList[[2]], "gp", list(cex=.8, col="blue"))
viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)

## ----viewTracksYlab,fig.cap='plot data with adjusted color and size of y-axis label',fig.width=8,fig.height=3----
setTrackStyleParam(trackList[[1]], "ylabgp", list(cex=.8, col="green"))
## set cex to avoid automatic adjust
setTrackStyleParam(trackList[[2]], "ylabgp", list(cex=.8, col="blue"))
setTrackStyleParam(trackList[[2]], "marginBottom", .2)
viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)

## ----viewTracksYlabTopBottom,fig.cap='plot data with adjusted y-axis label position',fig.width=8,fig.height=3----
setTrackStyleParam(trackList[[1]], "ylabpos", "bottomleft")
setTrackStyleParam(trackList[[2]], "ylabpos", "topright")
setTrackStyleParam(trackList[[2]], "marginTop", .2)
viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)

## ----viewTracksYlabUpsDown,fig.cap='plot data with adjusted transcripts name position',fig.width=8,fig.height=3----
trackListN <- trackList
setTrackStyleParam(trackListN[[3]], "ylabpos", "upstream")
setTrackStyleParam(trackListN[[4]], "ylabpos", "downstream")
## set cex to avoid automatic adjust
setTrackStyleParam(trackListN[[3]], "ylabgp", list(cex=.6))
setTrackStyleParam(trackListN[[4]], "ylabgp", list(cex=.6))
gr1 <- range(unname(unlist(GRangesList(sapply(trs, function(.ele) .ele$dat)))))
start(gr1) <- start(gr1) - 2000
end(gr1) <- end(gr1) + 2000
viewTracks(trackListN, gr=gr1, viewerStyle=viewerStyle)

## ----viewTracksCol,fig.cap='plot data with adjusted track color',fig.width=8,fig.height=3----
setTrackStyleParam(trackList[[1]], "color", c("green", "black"))
setTrackStyleParam(trackList[[2]], "color", c("black", "blue"))
for(i in 3:length(trackList)) 
    setTrackStyleParam(trackList[[i]], "color", "black")
viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)

## ----viewTracksHeight,fig.cap='plot data with adjusted track height',fig.width=8,fig.height=3----
trackListH <- trackList
setTrackStyleParam(trackListH[[1]], "height", .1)
setTrackStyleParam(trackListH[[2]], "height", .44)
for(i in 3:length(trackListH)){
    setTrackStyleParam(trackListH[[i]], "height", 
                       (1-(0.1+0.44))/(length(trackListH)-2))
}
viewTracks(trackListH, gr=gr, viewerStyle=viewerStyle)

## ----viewTracksNames,fig.cap='change the track names',fig.width=8,fig.height=3----
names(trackList) <- c("cpsf160", "fox2", rep("HSPA8", 5))
viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)

## ----viewTracksPaired,fig.cap='show two data in the same track',fig.width=8,fig.height=2.4----
cpsf160 <- importScore(file.path(extdata, "cpsf160.repA_-.wig"),
                       file.path(extdata, "cpsf160.repB_-.wig"),
                       format="WIG")
strand(cpsf160$dat) <- strand(cpsf160$dat2) <- "-"
setTrackStyleParam(cpsf160, "color", c("black", "red"))
viewTracks(trackList(trs, cpsf160), gr=gr, viewerStyle=viewerStyle)

## ----viewTracksFlipped,fig.cap='show data in the flipped track',fig.width=8,fig.height=4----
viewerStyleF <- viewerStyle
setTrackViewerStyleParam(viewerStyleF, "flip", TRUE)
setTrackViewerStyleParam(viewerStyleF, "xaxis", TRUE)
setTrackViewerStyleParam(viewerStyleF, "margin", c(.1, .05, .01, .01))
vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyleF)
addGuideLine(c(122929767, 122929969), vp=vp)
addArrowMark(list(x=122929650,
                  y=2),
             label="label",
             col="blue",
             vp=vp)

## ----themeBW,fig.cap='balck & white theme',fig.width=8,fig.height=4-----------
optSty <- optimizeStyle(trackList(repA, fox2, trs), theme="bw")
trackList <- optSty$tracks
viewerStyle <- optSty$style
vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)

## ----themeCol,fig.cap='colorful theme',fig.width=8,fig.height=4---------------
optSty <- optimizeStyle(trackList(repA, fox2, trs), theme="col")
trackList <- optSty$tracks
viewerStyle <- optSty$style
vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)

## ----themeSafe,fig.cap='safe theme',fig.width=8,fig.height=4------------------
optSty <- optimizeStyle(trackList(repA, fox2, trs), theme="safe")
trackList <- optSty$tracks
viewerStyle <- optSty$style
vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)

## ----axisBreak,fig.cap='axis with breaks',fig.width=8,fig.height=4------------
gr.breaks <- GRanges("chr11", 
                     IRanges(c(122929275, 122929575, 122929775), 
                             c(122929555, 122929725, 122930122)), 
                     strand="-", percentage=c(.4, .2, .4))
vp <- viewTracks(trackList, gr=gr.breaks, viewerStyle=viewerStyle)

## ----sessionInfo, results='asis'----------------------------------------------
sessionInfo()

