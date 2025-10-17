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

## ----out.width="50%", fig.asp=.6----------------------------------------------
plotGeneCount(curve = crv, clusters = celltype, 
              title = "Colored by cell type")

## ----eval=FALSE---------------------------------------------------------------
#  ### Based on Slingshot object
#  set.seed(6)
#  icMat <- evaluateK(counts = counts, sds = crv, k = 3:7, nGenes = 100,
#                     verbose = FALSE, plot = TRUE)
#  print(icMat[1:2, ])
#  
#  ### Downstream of any trajectory inference method using pseudotime and cell weights
#  set.seed(7)
#  pseudotime <- slingPseudotime(crv, na=FALSE)
#  cellWeights <- slingCurveWeights(crv)
#  icMat2 <- evaluateK(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
#                     k=3:7, nGenes = 100, verbose = FALSE, plot = TRUE)

## -----------------------------------------------------------------------------
### Based on Slingshot object
set.seed(6)
sce <- fitGAM(counts = counts, sds = crv, nknots = 6, verbose = FALSE)

### Downstream of any trajectory inference method using pseudotime and cell weights
set.seed(7)
pseudotime <- slingPseudotime(crv, na = FALSE)
cellWeights <- slingCurveWeights(crv)
sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
                 nknots = 6, verbose = FALSE)

## -----------------------------------------------------------------------------
batch <- factor(rep(c("A", "B"), each = ncol(counts)/2))
U <- model.matrix(~batch)
sceBatch <- fitGAM(counts = counts,
                   U = U,
                   sds = crv, 
                   nknots = 6, 
                   verbose = FALSE)


## -----------------------------------------------------------------------------
BPPARAM <- BiocParallel::bpparam()
BPPARAM # lists current options
BPPARAM$workers <- 2 # use 2 cores

sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
                 nknots = 6, verbose = FALSE, parallel=TRUE, BPPARAM = BPPARAM)

## -----------------------------------------------------------------------------
sce25 <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
                 nknots = 6, verbose = FALSE, genes = 1:25)

## -----------------------------------------------------------------------------
library(mgcv)
control <- gam.control()
control$maxit <- 1000 #set maximum number of iterations to 1K
# pass to control argument of fitGAM as below:
# 
# gamList <- fitGAM(counts = counts,
#                   pseudotime = slingPseudotime(crv, na = FALSE),
#                   cellWeights = slingCurveWeights(crv),
#                   control = control)

## -----------------------------------------------------------------------------
gamList <- fitGAM(counts,
                  pseudotime = slingPseudotime(crv, na = FALSE),
                  cellWeights = slingCurveWeights(crv),
                  nknots = 6, sce = FALSE)

## -----------------------------------------------------------------------------
summary(gamList[["Irf8"]])

## ----eval=FALSE---------------------------------------------------------------
#  pvalLineage <- getSmootherPvalues(gamList)
#  statLineage <- getSmootherTestStats(gamList)

## -----------------------------------------------------------------------------
sessionInfo()

