## ----echo = FALSE-------------------------------------------------------------
library(knitr)

## ----eval = FALSE-------------------------------------------------------------
#  if(!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
#  }
#  BiocManager::install("tradeSeq")

## ----warning=F, message=F-----------------------------------------------------
library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)

# For reproducibility
RNGversion("3.5.0")
palette(brewer.pal(8, "Dark2"))
data(countMatrix, package = "tradeSeq")
counts <- as.matrix(countMatrix)
rm(countMatrix)
data(crv, package = "tradeSeq")
data(celltype, package = "tradeSeq")

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(5)
#  icMat <- evaluateK(counts = counts, sds = crv, k = 3:10,
#                     nGenes = 200, verbose = T)

## ----echo=FALSE---------------------------------------------------------------
knitr::include_graphics("evalK_paul.png")

## -----------------------------------------------------------------------------
set.seed(7)
pseudotime <- slingPseudotime(crv, na = FALSE)
cellWeights <- slingCurveWeights(crv)
sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
                 nknots = 6, verbose = FALSE)

## -----------------------------------------------------------------------------
table(rowData(sce)$tradeSeq$converged)

## -----------------------------------------------------------------------------
assoRes <- associationTest(sce)
head(assoRes)

## -----------------------------------------------------------------------------
startRes <- startVsEndTest(sce)

## ----out.width="40%", fig.asp=1-----------------------------------------------
oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[3]]
plotSmoothers(sce, counts, gene = sigGeneStart)

## ----out.width="50%", fig.asp=.5----------------------------------------------
plotGeneCount(crv, counts, gene = sigGeneStart)

## -----------------------------------------------------------------------------
customRes <- startVsEndTest(sce, pseudotimeValues = c(0.1, 0.8))

## -----------------------------------------------------------------------------
endRes <- diffEndTest(sce)

## ----out.width="40%", fig.asp=1-----------------------------------------------
o <- order(endRes$waldStat, decreasing = TRUE)
sigGene <- names(sce)[o[1]]
plotSmoothers(sce, counts, sigGene)

## ----out.width="50%", fig.asp=.5----------------------------------------------
plotGeneCount(crv, counts, gene = sigGene)

## ----out.width="40%", fig.asp=1-----------------------------------------------
patternRes <- patternTest(sce)
oPat <- order(patternRes$waldStat, decreasing = TRUE)
head(rownames(patternRes)[oPat])
plotSmoothers(sce, counts, gene = rownames(patternRes)[oPat][4])

## ----out.width="50%", fig.asp=.5----------------------------------------------
plotGeneCount(crv, counts, gene = rownames(patternRes)[oPat][4])

## ----out.width="50%", fig.asp=.8----------------------------------------------
library(ggplot2)

patternRes$Gene <- rownames(patternRes)
patternRes$pattern <- patternRes$waldStat
patternRes <- patternRes[, c("Gene", "pattern")]

endRes$Gene <- rownames(endRes)
endRes$end <- endRes$waldStat
endRes <- endRes[, c("Gene", "end")]

compare <- merge(patternRes, endRes, by = "Gene", all = FALSE)
compare$transientScore <- 
  rank(-compare$end, ties.method = "min")^2 + rank(compare$pattern, ties.method = "random")^2

ggplot(compare, aes(x = log(pattern), y = log(end))) +
  geom_point(aes(col = transientScore)) +
  labs(x = "patternTest Wald Statistic (log scale)",
       y = "diffEndTest Wald Statistic (log scale)") +
  scale_color_continuous(low = "yellow", high = "red") +
  theme_classic()

## ----out.width="40%", fig.asp=1-----------------------------------------------
topTransient <- compare[which.max(compare$transientScore), "Gene"]
plotSmoothers(sce, counts, gene = topTransient)

## ----out.width="50%", fig.asp=.5----------------------------------------------
plotGeneCount(crv, counts, gene = topTransient)

## -----------------------------------------------------------------------------
head(
  compare[order(compare$transientScore, decreasing = TRUE), "Gene"],
  n = 5
)

## ----out.width="40%", fig.asp=1-----------------------------------------------
plotSmoothers(sce, counts, gene = "Irf8")

## ----out.width="50%", fig.asp=.5----------------------------------------------
plotGeneCount(crv, counts, gene = "Irf8")

## ----out.width="50%", fig.asp=.5----------------------------------------------
plotGeneCount(curve = crv, counts = counts,
              clusters = apply(slingClusterLabels(crv), 1, which.max),
              models = sce)


earlyDERes <- earlyDETest(sce, knots = c(1, 2))
oEarly <- order(earlyDERes$waldStat, decreasing = TRUE)
head(rownames(earlyDERes)[oEarly])

## ----out.width="40%", fig.asp=1-----------------------------------------------
plotSmoothers(sce, counts, gene = rownames(earlyDERes)[oEarly][2])

## ----out.width="50%", fig.asp=.5----------------------------------------------
plotGeneCount(crv, counts, gene = rownames(earlyDERes)[oEarly][2])

## -----------------------------------------------------------------------------
# testing against fold change threshold of 2
start2 <- startVsEndTest(sce, l2fc = log2(2))
# testing against fold change threshold of 1.5
pat2 <- patternTest(sce, l2fc = log2(1.5))

## -----------------------------------------------------------------------------
yhat <- predictCells(models = sce, gene = "Irf8")

ysmooth <- predictSmooth(models = sce, gene = "Irf8", nPoints = 40)

## ----clustering, warning=FALSE,message=F--------------------------------------
library(clusterExperiment)
nPointsClus <- 20
clusPat <- clusterExpressionPatterns(sce, nPoints = nPointsClus,
                                     genes = rownames(counts)[1:100])
clusterLabels <- primaryCluster(clusPat$rsec)

## ----plot clustering, out.width="25%", fig.asp=1, fig.show="hold"-------------
cUniq <- unique(clusterLabels)
cUniq <- cUniq[!cUniq == -1] # remove unclustered genes

for (xx in cUniq[1:4]) {
  cId <- which(clusterLabels == xx)
  p <- ggplot(data = data.frame(x = 1:nPointsClus,
                                y = rep(range(clusPat$yhatScaled[cId, ]),
                                        nPointsClus / 2)),
              aes(x = x, y = y)) +
    geom_point(alpha = 0) +
    labs(title = paste0("Cluster ", xx),  x = "Pseudotime", y = "Normalized expression") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  for (ii in 1:length(cId)) {
    geneId <- rownames(clusPat$yhatScaled)[cId[ii]]
    p <- p +
      geom_line(data = data.frame(x = rep(1:nPointsClus, 2),
                                  y = clusPat$yhatScaled[geneId, ],
                                  lineage = rep(0:1, each = nPointsClus)),
                aes(col = as.character(lineage), group = lineage), lwd = 1.5)
  }
  p <- p + guides(color = FALSE) +
    scale_color_manual(values = c("orange", "darkseagreen3"),
                       breaks = c("0", "1"))  
  print(p)
}

## -----------------------------------------------------------------------------
sessionInfo()

