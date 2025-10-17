## ----echo=FALSE, results="hide", warning=FALSE--------------------------------
suppressPackageStartupMessages({
  library(trackViewer)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(InteractionSet)
})
knitr::opts_chunk$set(warning=FALSE, message=FALSE)

## ----plotback2back------------------------------------------------------------
library(trackViewer)
library(InteractionSet)
gi <- readRDS(system.file("extdata", "nij.chr6.51120000.53200000.gi.rds", package="trackViewer"))
head(gi)
## hicexplorer:hicConvertFormat tool can be used to convert other formats into GInteractions
## eg: hicConvertFormat -m mESC_rep.hic --inputFormat hic --outputFormat cool -o mESC_rep.mcool
##     hicConvertFormat -m mESC_rep.mcool::resolutions/10000 --inputFormat cool --outputFormat ginteractions -o mESC_rep.ginteractions --resolutions 10000
## please note that metadata:score is used for plot.
gi$border_color <- NA ## highlight some regions
gi$border_color[sample(seq_along(gi), 20)] <- sample(1:7, 20, replace=TRUE)
## The TADs will be drawn as lines at points start(first), center point, end(second).
tads <- GInteractions(
  GRanges("chr6", 
          IRanges(c(51130001, 51130001, 51450001, 52210001), width = 20000)),
  GRanges("chr6", 
          IRanges(c(51530001, 52170001, 52210001, 53210001), width = 20000)))
range <- GRanges("chr6", IRanges(51120000, 53200000))
heatmap <- gi2track(gi)
ctcf <- readRDS(system.file("extdata", "ctcf.sample.rds", package="trackViewer"))
viewTracks(trackList(ctcf, heatmap, heightDist = c(1, 3)), 
           gr=range, autoOptimizeStyle = TRUE)
## add TAD information
addInteractionAnnotation(tads, "heatmap", grid.lines, gp=gpar(col="#E69F00", lwd=3, lty=3))
## add highlight interested regions
gi_sub <- gi[order(gi$score, decreasing = TRUE)]
gi_sub <- head(gi_sub[distance(first(gi_sub), second(gi_sub))>200000], n=5)
start(regions(gi_sub)) <- start(regions(gi_sub))-40000
end(regions(gi_sub)) <- end(regions(gi_sub))+40000
addInteractionAnnotation(gi_sub, "heatmap", grid.polygon, gp=gpar(col="red", lwd=2, lty=2, fill=NA))
## add interesting anchor at giving coordinate.
addInteractionAnnotation(52900000, "heatmap", gp=gpar(col="blue", lwd=3))
addInteractionAnnotation(-52900000, "heatmap", gp=gpar(col="cyan", lwd=3, lty=4))
## view the interaction data back to back.
## Please make sure the data are normalized.
gi2 <- gi
set.seed(123)
gi2$score <- gi$score + rnorm(length(gi), sd = sd(gi$score))
back2back <- gi2track(gi, gi2)
## change the color
setTrackStyleParam(back2back, "breaks", 
                   c(seq(from=0, to=50, by=10), 200))
setTrackStyleParam(back2back, "color",
                   c("lightblue", "yellow", "red"))
## chang the lim of y-axis (by default, [0, 1])
setTrackStyleParam(back2back, "ylim", c(0, .5))
viewTracks(trackList(ctcf, back2back, heightDist=c(1, 5)),
           gr=range, autoOptimizeStyle = TRUE)
addInteractionAnnotation(tads, "back2back", grid.lines,
                         gp=gpar(col="cyan", lwd=3, lty=2))
addInteractionAnnotation(-52208000, "back2back", gp=gpar(col="blue", lwd=3),
                         panel="top")
addInteractionAnnotation(51508000, "back2back", gp=gpar(col="gray", lwd=3, lty=2),
                         panel="bottom")

## ----plotLinks, fig.width=6, fig.height=3-------------------------------------
setTrackStyleParam(heatmap, "tracktype", "link")
setTrackStyleParam(heatmap, "breaks", 
                   c(seq(from=0, to=50, by=10), 200))
setTrackStyleParam(heatmap, "color",
                   c("lightblue", "yellow", "red"))
## filter the links to simulate the real data
keep <- distance(heatmap$dat, heatmap$dat2) > 5e5 & heatmap$dat$score>20
heatmap$dat <- heatmap$dat[keep]
heatmap$dat2 <- heatmap$dat2[keep]
viewTracks(trackList(heatmap), gr=range, autoOptimizeStyle = TRUE)

## ----inportHic----------------------------------------------------------------
hic <- system.file("extdata", "test_chr22.hic", package = "trackViewer",
                    mustWork=TRUE)
if(.Platform$OS.type!="windows"){
importGInteractions(file=hic, format="hic",
                    ranges=GRanges("22", IRanges(50000000, 100000000)),
                    out = "GInteractions")
}

## ----importCool, eval=FALSE---------------------------------------------------
#  cool <- system.file("extdata", "test.mcool", package = "trackViewer",
#                       mustWork=TRUE)
#  importGInteractions(file=cool, format="cool",
#                      resolution = 2,
#                      ranges=GRanges("chr1", IRanges(10, 28)),
#                      out = "GInteractions")

## ----plotGInteractions--------------------------------------------------------
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(InteractionSet)
gi <- readRDS(system.file("extdata", "gi.rds", package="trackViewer"))
range <- GRanges("chr2", IRanges(234500000, 235000000))
feature.gr <- suppressMessages(genes(TxDb.Hsapiens.UCSC.hg19.knownGene))
feature.gr <- subsetByOverlaps(feature.gr, regions(gi))
feature.gr$col <- sample(1:7, length(feature.gr), replace=TRUE)
feature.gr$type <- sample(c("promoter", "enhancer", "gene"), 
                          length(feature.gr), replace=TRUE, 
                          prob=c(0.1, 0.2, 0.7))
plotGInteractions(gi, range, feature.gr)

## ----sessionInfo, results='asis'----------------------------------------------
sessionInfo()

