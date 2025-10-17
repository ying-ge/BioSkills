## ----options, results="hide", include=FALSE, cache=FALSE, message=FALSE-------
knitr::opts_chunk$set(fig.align="center", cache=TRUE,error=FALSE, #stop on error
fig.width=5, fig.height=5, autodep=TRUE,
results="markup", echo=TRUE, eval=TRUE)
#knitr::opts_knit$set(stop_on_error = 2L) #really make it stop
#knitr::dep_auto()
options(getClass.msg=FALSE)
graphics:::par(pch = 16, las = 1)
set.seed(12345) ## for reproducibility
library(SingleCellExperiment)
library(slingshot, quietly = TRUE)

## ----dataSetup_sim------------------------------------------------------------
# generate synthetic count data representing a single lineage
means <- rbind(
    # non-DE genes
    matrix(rep(rep(c(0.1,0.5,1,2,3), each = 300),100),
        ncol = 300, byrow = TRUE),
    # early deactivation
    matrix(rep(exp(atan( ((300:1)-200)/50 )),50), ncol = 300, byrow = TRUE),
    # late deactivation
    matrix(rep(exp(atan( ((300:1)-100)/50 )),50), ncol = 300, byrow = TRUE),
    # early activation
    matrix(rep(exp(atan( ((1:300)-100)/50 )),50), ncol = 300, byrow = TRUE),
    # late activation
    matrix(rep(exp(atan( ((1:300)-200)/50 )),50), ncol = 300, byrow = TRUE),
    # transient
    matrix(rep(exp(atan( c((1:100)/33, rep(3,100), (100:1)/33) )),50), 
        ncol = 300, byrow = TRUE)
)
counts <- apply(means,2,function(cell_means){
    total <- rnbinom(1, mu = 7500, size = 4)
    rmultinom(1, total, cell_means)
})
rownames(counts) <- paste0('G',1:750)
colnames(counts) <- paste0('c',1:300)
sce <- SingleCellExperiment(assays = List(counts = counts))

## ----data_sling---------------------------------------------------------------
library(slingshot, quietly = FALSE)
data("slingshotExample")
rd <- slingshotExample$rd
cl <- slingshotExample$cl

dim(rd) # data representing cells in a reduced dimensional space
length(cl) # vector of cluster labels

## ----genefilt-----------------------------------------------------------------
# filter genes down to potential cell-type markers
# at least M (15) reads in at least N (15) cells
geneFilter <- apply(assays(sce)$counts,1,function(x){
    sum(x >= 3) >= 10
})
sce <- sce[geneFilter, ]

## ----norm---------------------------------------------------------------------
FQnorm <- function(counts){
    rk <- apply(counts,2,rank,ties.method='min')
    counts.sort <- apply(counts,2,sort)
    refdist <- apply(counts.sort,1,median)
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
}
assays(sce)$norm <- FQnorm(assays(sce)$counts)

## ----pca, cache=TRUE----------------------------------------------------------
pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]

plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)

## ----umap, cache=TRUE---------------------------------------------------------
library(uwot)
rd2 <- uwot::umap(t(log1p(assays(sce)$norm)))
colnames(rd2) <- c('UMAP1', 'UMAP2')

plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)

## ----add_RDs, cache=TRUE------------------------------------------------------
reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)

## ----clustering_mclust--------------------------------------------------------
library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification
colData(sce)$GMM <- cl1

library(RColorBrewer)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

## ----clustering---------------------------------------------------------------
cl2 <- kmeans(rd1, centers = 4)$cluster
colData(sce)$kmeans <- cl2

plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)

## ----sling_sce----------------------------------------------------------------
sce <- slingshot(sce, clusterLabels = 'GMM', reducedDim = 'PCA')

## ----plot_curve_1-------------------------------------------------------------
summary(sce$slingPseudotime_1)

library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

## ----plot_curve_2-------------------------------------------------------------
plot(reducedDims(sce)$PCA, col = brewer.pal(9,'Set1')[sce$GMM], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')

## ----tradeseq, eval=FALSE-----------------------------------------------------
#  library(tradeSeq)
#  
#  # fit negative binomial GAM
#  sce <- fitGAM(sce)
#  
#  # test for dynamic expression
#  ATres <- associationTest(sce)

## ----heatmaps, eval=FALSE-----------------------------------------------------
#  topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:250]
#  pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
#  heatdata <- assays(sce)$counts[topgenes, pst.ord]
#  heatclus <- sce$GMM[pst.ord]
#  
#  heatmap(log1p(heatdata), Colv = NA,
#          ColSideColors = brewer.pal(9,"Set1")[heatclus])

## ----heatmapsREAL, echo=FALSE, fig.height=7-----------------------------------
topgenes <- paste0('G',501:750)
# tradeSeq has too many dependencies (174 at the time of this writing), but I 
# promise I actually ran it and got this result. This is a *very* clean example 
# dataset
pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
heatdata <- assays(sce)$counts[topgenes, pst.ord]
heatclus <- sce$GMM[pst.ord]
heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])

## ----sling_lines_unsup--------------------------------------------------------
lin1 <- getLineages(rd, cl, start.clus = '1')
lin1
plot(rd, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
lines(SlingshotDataSet(lin1), lwd = 3, col = 'black')

## ----lines_sup_end------------------------------------------------------------
lin2 <- getLineages(rd, cl, start.clus= '1', end.clus = '3')

plot(rd, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
lines(SlingshotDataSet(lin2), lwd = 3, col = 'black', show.constraints = TRUE)

## ----curves-------------------------------------------------------------------
crv1 <- getCurves(lin1)
crv1
plot(rd, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
lines(SlingshotDataSet(crv1), lwd = 3, col = 'black')

## ----sling_approxpoints-------------------------------------------------------
sce5 <- slingshot(sce, clusterLabels = 'GMM', reducedDim = 'PCA',
                   approx_points = 5)

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce5$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce5)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce5), lwd=2, col='black')

## ----sling_omega--------------------------------------------------------------
rd2 <- rbind(rd, cbind(rd[,2]-12, rd[,1]-6))
cl2 <- c(cl, cl + 10)
pto2 <- slingshot(rd2, cl2, omega = TRUE, start.clus = c(1,11))

plot(rd2, pch=16, asp = 1,
     col = c(brewer.pal(9,"Set1"), brewer.pal(8,"Set2"))[cl2])
lines(SlingshotDataSet(pto2), type = 'l', lwd=2, col='black')

## ----sling_multtraj-----------------------------------------------------------
plot(rd2, pch=16, asp = 1,
     col = c(brewer.pal(9,"Set1"), brewer.pal(8,"Set2"))[cl2])
lines(SlingshotDataSet(pto2), lwd=2, col='black')

## ----sling_predict------------------------------------------------------------
# our original PseudotimeOrdering
pto <- sce$slingshot

# simulate new cells in PCA space
newPCA <- reducedDim(sce, 'PCA') + rnorm(2*ncol(sce), sd = 2)

# project onto trajectory
newPTO <- slingshot::predict(pto, newPCA)

## ----plot_sling_predict-------------------------------------------------------
newplotcol <- colors[cut(slingPseudotime(newPTO)[,1], breaks=100)]
plot(reducedDims(sce)$PCA, col = 'grey', bg = 'grey', pch=21, asp = 1,
     xlim = range(newPCA[,1]), ylim = range(newPCA[,2]))
lines(SlingshotDataSet(sce), lwd=2, col = 'black')
points(slingReducedDim(newPTO), col = newplotcol, pch = 16)

## ----session------------------------------------------------------------------
sessionInfo()

