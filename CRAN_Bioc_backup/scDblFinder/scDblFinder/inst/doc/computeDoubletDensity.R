## ----echo=FALSE, message=FALSE------------------------------------------------
knitr::opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)
library(BiocStyle)

## -----------------------------------------------------------------------------
library(scRNAseq)
sce <- BachMammaryData(samples="G_1")

set.seed(1001)
sce <- sce[,sample(ncol(sce), 1000)]

## -----------------------------------------------------------------------------
library(scuttle)
sce <- logNormCounts(sce)

library(scran)
dec <- modelGeneVar(sce)
hvgs <- getTopHVGs(dec, n=1000)

library(scater)
set.seed(1002)
sce <- runPCA(sce, ncomponents=10, subset_row=hvgs)
sce <- runTSNE(sce, dimred="PCA")

## -----------------------------------------------------------------------------
set.seed(1003)
library(scDblFinder)
scores <- computeDoubletDensity(sce, subset.row=hvgs)
plotTSNE(sce, colour_by=I(log1p(scores)))

## ----echo=FALSE---------------------------------------------------------------
# Sanity check that the plot has one cluster with much higher scores.
# If this fails, we probably need to pick a more demonstrative example.
library(bluster)
clusters <- clusterRows(reducedDim(sce, "PCA"), NNGraphParam())
by.clust <- split(scores, clusters)
med.scores <- sort(vapply(by.clust, median, 0), decreasing=TRUE)
stopifnot(med.scores[1] > med.scores[2] * 4)

## -----------------------------------------------------------------------------
sessionInfo()

