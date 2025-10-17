## ----echo=FALSE, results="hide", message=FALSE--------------------------------
knitr::opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)
library(BiocStyle)

## -----------------------------------------------------------------------------
library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()
hpca.se

## -----------------------------------------------------------------------------
library(scRNAseq)
hESCs <- LaMannoBrainData('human-es')
hESCs <- hESCs[,1:100]

## -----------------------------------------------------------------------------
library(SingleR)
pred.hesc <- SingleR(test = hESCs, ref = hpca.se, assay.type.test=1,
    labels = hpca.se$label.main)

## -----------------------------------------------------------------------------
pred.hesc
# Summarizing the distribution:
table(pred.hesc$labels)

## -----------------------------------------------------------------------------
library(scRNAseq)
sceM <- MuraroPancreasData()

# One should normally do cell-based quality control at this point, but for
# brevity's sake, we will just remove the unlabelled libraries here.
sceM <- sceM[,!is.na(sceM$label)]

# SingleR() expects reference datasets to be normalized and log-transformed.
library(scuttle)
sceM <- logNormCounts(sceM)

## -----------------------------------------------------------------------------
sceG <- GrunPancreasData()
sceG <- sceG[,colSums(counts(sceG)) > 0] # Remove libraries with no counts.
sceG <- logNormCounts(sceG) 

## -----------------------------------------------------------------------------
pred.grun <- SingleR(test=sceG, ref=sceM, labels=sceM$label, de.method="wilcox")
table(pred.grun$labels)

## -----------------------------------------------------------------------------
plotScoreHeatmap(pred.grun)

## -----------------------------------------------------------------------------
plotDeltaDistribution(pred.grun, ncol = 3)

## -----------------------------------------------------------------------------
summary(is.na(pred.grun$pruned.labels))

## -----------------------------------------------------------------------------
all.markers <- metadata(pred.grun)$de.genes
sceG$labels <- pred.grun$labels

# Beta cell-related markers
library(scater)
plotHeatmap(sceG, order_columns_by="labels",
    features=unique(unlist(all.markers$beta))) 

## -----------------------------------------------------------------------------
sessionInfo()

