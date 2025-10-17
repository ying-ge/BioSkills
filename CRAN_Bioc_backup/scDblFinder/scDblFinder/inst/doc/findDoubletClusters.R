## ----echo=FALSE, message=FALSE------------------------------------------------
knitr::opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)
library(BiocStyle)

## -----------------------------------------------------------------------------
library(scRNAseq)
sce <- BachMammaryData(samples="G_2")

set.seed(1000)
sce <- sce[,sample(ncol(sce), 500)]

## -----------------------------------------------------------------------------
library(scuttle)
sce <- logNormCounts(sce)

library(scran)
dec <- modelGeneVar(sce)

library(scater)
set.seed(1000)
sce <- runPCA(sce, ncomponents=10, subset_row=getTopHVGs(dec, n=1000))

library(bluster)
clusters <- clusterRows(reducedDim(sce, "PCA"), NNGraphParam())

sce <- runTSNE(sce, dimred="PCA")
plotTSNE(sce, colour_by=I(clusters), text_by=I(clusters))

## -----------------------------------------------------------------------------
library(scDblFinder)
tab <- findDoubletClusters(sce, clusters)
tab

## ----echo=FALSE---------------------------------------------------------------
# Sanity check that one of the clusters is a good doublet candidate.
# If this fails, we probably need to pick a more demonstrative example.
stopifnot(rownames(tab)[1]=="6")
stopifnot(tab[1,"num.de"]==0)

## -----------------------------------------------------------------------------
sessionInfo()

