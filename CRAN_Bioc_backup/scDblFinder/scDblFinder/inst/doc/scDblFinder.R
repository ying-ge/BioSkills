## ----include=FALSE------------------------------------------------------------
library(BiocStyle)

## ----eval=FALSE---------------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("scDblFinder")
#  
#  # or, to get that latest developments:
#  BiocManager::install("plger/scDblFinder")

## ----warning=FALSE------------------------------------------------------------
set.seed(123)
suppressPackageStartupMessages(library(scDblFinder))
# we create a dummy dataset; since it's small we set a higher doublet rate
sce <- mockDoubletSCE(dbl.rate=0.1, ngenes=300 )
# we run scDblFinder (providing the unusually high doublet rate)
sce <- scDblFinder(sce, dbr=0.1)

## -----------------------------------------------------------------------------
table(truth=sce$type, call=sce$scDblFinder.class)

## -----------------------------------------------------------------------------
sce <- scDblFinder(sce, clusters="cluster")
table(truth=sce$type, call=sce$scDblFinder.class)

## ----eval=FALSE---------------------------------------------------------------
#  library(BiocParallel)
#  sce <- scDblFinder(sce, samples="sample_id", BPPARAM=MulticoreParam(3))
#  table(sce$scDblFinder.class)

## -----------------------------------------------------------------------------
sessionInfo()

