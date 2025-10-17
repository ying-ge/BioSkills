## ----echo=FALSE, results="hide", warning=FALSE--------------------------------
suppressPackageStartupMessages({
    library(trackViewer)
    library(rtracklayer)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(org.Hs.eg.db)
    library(VariantAnnotation)
  library(httr)
})
knitr::opts_chunk$set(warning=FALSE, message=FALSE)

## ----fig.width=4.5,fig.height=3-----------------------------------------------
library(trackViewer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(rtracklayer)
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

## ----fig.width=4.5,fig.height=3-----------------------------------------------
methy$color <- 3
methy$border <- "gray"
## Score info is required and the score must be a number in [0, 1]
m <- max(methy$score)
methy$score <- methy$score/m
dandelion.plot(methy, features, ranges=gr, type="fan")

## ----fig.width=4.5,fig.height=4-----------------------------------------------
methy$color <- rep(list(c(3, 5)), length(methy))
methy$score2 <- (max(methy$score) - methy$score)/m
legends <- list(list(labels=c("s1", "s2"), fill=c(3, 5)))
dandelion.plot(methy, features, ranges=gr, type="pie", legend=legends)

## ----fig.width=4.5,fig.height=3.5---------------------------------------------
## Less dandelions
dandelion.plot(methy, features, ranges=gr, type="circle", maxgaps=1/10)

## ----fig.width=4.5,fig.height=4-----------------------------------------------
## More dandelions
dandelion.plot(methy, features, ranges=gr, type="circle", maxgaps=1/100)

## ----fig.width=4,fig.height=4-------------------------------------------------
maxgaps <- tile(gr, n = 10)[[1]]
dandelion.plot(methy, features, ranges=gr, type="circle", maxgaps=maxgaps)

## ----fig.width=4.5,fig.height=3-----------------------------------------------
dandelion.plot(methy, features, ranges=gr, type="pie", 
               maxgaps=1/100, yaxis = TRUE, heightMethod = mean,
               ylab='mean of methy scores')

## ----fig.width=4.5,fig.height=3-----------------------------------------------
yaxis = c(0, 0.5, 1)
dandelion.plot(methy, features, ranges=gr, type="pie", 
               maxgaps=1/100, yaxis = yaxis, heightMethod = mean,
               ylab='mean of methy scores')

## ----sessionInfo, results='asis'----------------------------------------------
sessionInfo()

