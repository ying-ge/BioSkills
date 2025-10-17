## ----echo=FALSE, results="hide", message=FALSE--------------------------------
require(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)

## ----setup, echo=FALSE, message=FALSE-----------------------------------------
library(scran)
library(BiocParallel)
register(SerialParam()) # avoid problems with fastMNN parallelization.
set.seed(100)

## -----------------------------------------------------------------------------
library(scRNAseq)
sce <- GrunPancreasData()
sce

## -----------------------------------------------------------------------------
library(scuttle)
qcstats <- perCellQCMetrics(sce)
qcfilter <- quickPerCellQC(qcstats, percent_subsets="altexps_ERCC_percent")
sce <- sce[,!qcfilter$discard]
summary(qcfilter$discard)

## -----------------------------------------------------------------------------
library(scran)
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters=clusters)
summary(sizeFactors(sce))
sce <- logNormCounts(sce)

## -----------------------------------------------------------------------------
dec <- modelGeneVar(sce)
plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
curve(metadata(dec)$trend(x), col="blue", add=TRUE)

## -----------------------------------------------------------------------------
dec2 <- modelGeneVarWithSpikes(sce, 'ERCC')
plot(dec2$mean, dec2$total, xlab="Mean log-expression", ylab="Variance")
points(metadata(dec2)$mean, metadata(dec2)$var, col="red")
curve(metadata(dec2)$trend(x), col="blue", add=TRUE)

## ----fig.wide=TRUE, fig.asp=1.5-----------------------------------------------
# Turned off weighting to avoid overfitting for each donor.
dec3 <- modelGeneVar(sce, block=sce$donor, density.weights=FALSE)
per.block <- dec3$per.block
par(mfrow=c(3, 2))
for (i in seq_along(per.block)) {
    decX <- per.block[[i]]
    plot(decX$mean, decX$total, xlab="Mean log-expression", 
        ylab="Variance", main=names(per.block)[i])
    curve(metadata(decX)$trend(x), col="blue", add=TRUE)
}

## -----------------------------------------------------------------------------
# Get the top 10% of genes.
top.hvgs <- getTopHVGs(dec, prop=0.1)

# Get the top 2000 genes.
top.hvgs2 <- getTopHVGs(dec, n=2000)

# Get all genes with positive biological components.
top.hvgs3 <- getTopHVGs(dec, var.threshold=0)

# Get all genes with FDR below 5%.
top.hvgs4 <- getTopHVGs(dec, fdr.threshold=0.05)

## -----------------------------------------------------------------------------
sce <- fixedPCA(sce, subset.row=top.hvgs)
reducedDimNames(sce)

## -----------------------------------------------------------------------------
sced <- denoisePCA(sce, dec2, subset.row=getTopHVGs(dec2, prop=0.1))
ncol(reducedDim(sced, "PCA"))

## -----------------------------------------------------------------------------
output <- getClusteredPCs(reducedDim(sce))
npcs <- metadata(output)$chosen
reducedDim(sce, "PCAsub") <- reducedDim(sce, "PCA")[,1:npcs,drop=FALSE]
npcs

## -----------------------------------------------------------------------------
# In this case, using the PCs that we chose from getClusteredPCs().
g <- buildSNNGraph(sce, use.dimred="PCAsub")
cluster <- igraph::cluster_walktrap(g)$membership

# Assigning to the 'colLabels' of the 'sce'.
colLabels(sce) <- factor(cluster)
table(colLabels(sce))

## -----------------------------------------------------------------------------
library(scater)
sce <- runTSNE(sce, dimred="PCAsub")
plotTSNE(sce, colour_by="label", text_by="label")

## -----------------------------------------------------------------------------
library(bluster)
ratio <- pairwiseModularity(g, cluster, as.ratio=TRUE)

library(pheatmap)
pheatmap(log10(ratio+1), cluster_cols=FALSE, cluster_rows=FALSE,
    col=rev(heat.colors(100)))

## -----------------------------------------------------------------------------
ass.prob <- bootstrapStability(sce, FUN=function(x) {
    g <- buildSNNGraph(x, use.dimred="PCAsub")
    igraph::cluster_walktrap(g)$membership
}, clusters=sce$cluster)

pheatmap(ass.prob, cluster_cols=FALSE, cluster_rows=FALSE,
    col=colorRampPalette(c("white", "blue"))(100))

## -----------------------------------------------------------------------------
subout <- quickSubCluster(sce, groups=colLabels(sce))
table(metadata(subout)$subcluster) # formatted as '<parent>.<subcluster>'

## -----------------------------------------------------------------------------
# Uses clustering information from 'colLabels(sce)' by default:
markers <- scoreMarkers(sce)
markers
colnames(markers[[1]])

## -----------------------------------------------------------------------------
# Just showing the first few columns for brevity.
markers[[1]][order(markers[[1]]$mean.AUC, decreasing=TRUE),1:4]

## -----------------------------------------------------------------------------
markers <- scoreMarkers(sce, full.stats=TRUE)
markers[[1]]$full.logFC.cohen

## -----------------------------------------------------------------------------
# Using the first 200 HVs, which are the most interesting anyway.
of.interest <- top.hvgs[1:200]
cor.pairs <- correlatePairs(sce, subset.row=of.interest)
cor.pairs

## -----------------------------------------------------------------------------
cor.pairs2 <- correlatePairs(sce, subset.row=of.interest, block=sce$donor)

## -----------------------------------------------------------------------------
cor.genes <- correlateGenes(cor.pairs)
cor.genes

## -----------------------------------------------------------------------------
y <- convertTo(sce, type="edgeR")

## -----------------------------------------------------------------------------
sessionInfo()

