## ----echo=FALSE, results="hide", message=FALSE--------------------------------
require(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)

## ----setup, echo=FALSE, message=FALSE-----------------------------------------
library(batchelor)

## -----------------------------------------------------------------------------
library(scRNAseq)
sce1 <- ZeiselBrainData()
sce1
sce2 <- TasicBrainData()
sce2

## -----------------------------------------------------------------------------
library(scuttle)
sce1 <- addPerCellQC(sce1, subsets=list(Mito=grep("mt-", rownames(sce1))))
qc1 <- quickPerCellQC(colData(sce1), sub.fields="subsets_Mito_percent")
sce1 <- sce1[,!qc1$discard]

sce2 <- addPerCellQC(sce2, subsets=list(Mito=grep("mt_", rownames(sce2))))
qc2 <- quickPerCellQC(colData(sce2), sub.fields="subsets_Mito_percent")
sce2 <- sce2[,!qc2$discard]

## -----------------------------------------------------------------------------
universe <- intersect(rownames(sce1), rownames(sce2))
sce1 <- sce1[universe,]
sce2 <- sce2[universe,]

## -----------------------------------------------------------------------------
out <- multiBatchNorm(sce1, sce2)
sce1 <- out[[1]]
sce2 <- out[[2]]

## -----------------------------------------------------------------------------
library(scran)
dec1 <- modelGeneVar(sce1)
dec2 <- modelGeneVar(sce2)
combined.dec <- combineVar(dec1, dec2)
chosen.hvgs <- getTopHVGs(combined.dec, n=5000)

## -----------------------------------------------------------------------------
combined <- correctExperiments(A=sce1, B=sce2, PARAM=NoCorrectParam())

library(scater)
set.seed(100)
combined <- runPCA(combined, subset_row=chosen.hvgs)
combined <- runTSNE(combined, dimred="PCA")
plotTSNE(combined, colour_by="batch")

## -----------------------------------------------------------------------------
library(batchelor)
set.seed(101)
f.out <- fastMNN(A=sce1, B=sce2, subset.row=chosen.hvgs)
str(reducedDim(f.out, "corrected"))

## -----------------------------------------------------------------------------
rle(f.out$batch)

## -----------------------------------------------------------------------------
set.seed(101)
f.out2 <- fastMNN(combined, batch=combined$batch, subset.row=chosen.hvgs)
str(reducedDim(f.out2, "corrected"))

## -----------------------------------------------------------------------------
set.seed(103)
f.out <- runTSNE(f.out, dimred="corrected")
plotTSNE(f.out, colour_by="batch")

## -----------------------------------------------------------------------------
cor.exp <- assay(f.out)[1,]
hist(cor.exp, xlab="Corrected expression for gene 1", col="grey80") 

## -----------------------------------------------------------------------------
# Using fewer genes as it is much, much slower. 
fewer.hvgs <- head(order(combined.dec$bio, decreasing=TRUE), 100)
classic.out <- mnnCorrect(sce1, sce2, subset.row=fewer.hvgs)

## -----------------------------------------------------------------------------
classic.out

## -----------------------------------------------------------------------------
# Removing the 'unclassified' cluster, which makes no sense:
not.unclass <- sce2$broad_type!="Unclassified"
clust.out <- clusterMNN(sce1, sce2[,not.unclass],
    subset.row=chosen.hvgs,
    clusters=list(sce1$level1class, sce2$broad_type[not.unclass])) 

## -----------------------------------------------------------------------------
clust.info <- metadata(clust.out)$cluster
split(clust.info$cluster, clust.info$meta)

## -----------------------------------------------------------------------------
rescale.out <- rescaleBatches(sce1, sce2)
rescale.out

## -----------------------------------------------------------------------------
rescale.out <- runPCA(rescale.out, subset_row=chosen.hvgs,
    exprs_values="corrected")
plotPCA(rescale.out, colour_by="batch")

## -----------------------------------------------------------------------------
regress.out <- regressBatches(sce1, sce2)
assay(regress.out)

## -----------------------------------------------------------------------------
# Pretend the first X cells in each batch are controls.
restrict <- list(1:100, 1:200) 
rescale.out <- rescaleBatches(sce1, sce2, restrict=restrict)

## -----------------------------------------------------------------------------
normed <- multiBatchNorm(A=sce1, B=sce2,
    norm.args=list(use_altexps=FALSE))
names(normed)

## -----------------------------------------------------------------------------
set.seed(100)

# Using the same BSPARAM argument as fastMNN(), for speed.
pca.out <- multiBatchPCA(A=sce1, B=sce2, subset.row=chosen.hvgs,
    BSPARAM=BiocSingular::IrlbaParam(deferred=TRUE))
names(pca.out)

## -----------------------------------------------------------------------------
sessionInfo()

