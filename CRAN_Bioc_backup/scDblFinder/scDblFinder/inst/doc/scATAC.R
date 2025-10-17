## ----include=FALSE------------------------------------------------------------
library(BiocStyle)

## -----------------------------------------------------------------------------
suppressPackageStartupMessages(library(scDblFinder))
# we use a dummy SingleCellExperiment as example:
sce <- mockDoubletSCE(ngenes=300)
# setting low number of artificial doublets (same as ncells) just for speedup:
sce <- scDblFinder(sce, artificialDoublets=1, aggregateFeatures=TRUE, nfeatures=25, processing="normFeatures")

## -----------------------------------------------------------------------------
# here we use a dummy fragment file for example:
fragfile <- system.file("extdata", "example_fragments.tsv.gz", package="scDblFinder")

# we might also give a GRanges of repeat elements, so that these regions are excluded:
suppressPackageStartupMessages(library(GenomicRanges))
repeats <- GRanges("chr6", IRanges(1000,2000))
# it's better to combine these with mitochondrial and sex chromosomes
otherChroms <- GRanges(c("M","chrM","MT","X","Y","chrX","chrY"),IRanges(1L,width=10^8))
# here since I don't know what chromosome notation you'll be using I've just put them all,
# although this will trigger a warning when combining them:
toExclude <- suppressWarnings(c(repeats, otherChroms))
# we then launch the method
res <- amulet(fragfile, regionsToExclude=toExclude)
res

## ----eval=FALSE---------------------------------------------------------------
#  # not run
#  d <- clamulet("path/to/fragments.tsv.gz")

## -----------------------------------------------------------------------------
d <- clamulet(fragfile, k=2, nfeatures=3)
d

## ----eval=FALSE---------------------------------------------------------------
#  res$scDblFinder.p <- 1-colData(sce)[row.names(res), "scDblFinder.score"]
#  res$combined <- apply(res[,c("scDblFinder.p", "p.value")], 1, FUN=function(x){
#    x[x<0.001] <- 0.001 # prevent too much skew from very small or 0 p-values
#    suppressWarnings(aggregation::fisher(x))
#  })

## -----------------------------------------------------------------------------
sessionInfo()

