## ----eval = FALSE, echo = TRUE, message = FALSE, warning = FALSE--------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("scds", version = "3.9")

## ----eval = FALSE, echo = TRUE, message = FALSE, warning = FALSE--------------
#  library(devtools)
#  devtools::install_github('kostkalab/scds')
#  

## ----prelims, eval = TRUE, echo = TRUE, message = FALSE, warning = FALSE------
library(scds)
library(scater)
library(rsvd)
library(Rtsne)
library(cowplot)
set.seed(30519)
data("sce_chcl")
sce = sce_chcl #- less typing
dim(sce)


## ----doublets, eval = TRUE, echo = TRUE, message = FALSE, warning = FALSE-----
table(sce$hto_classification_global)

## ----proj, eval = TRUE, echo = TRUE, message = FALSE, warning = FALSE---------
logcounts(sce) = log1p(counts(sce))
vrs            = apply(logcounts(sce),1,var)
pc             = rpca(t(logcounts(sce)[order(vrs,decreasing=TRUE)[1:100],]))
ts             = Rtsne(pc$x[,1:10],verb=FALSE)

reducedDim(sce,"tsne") = ts$Y; rm(ts,vrs,pc)
plotReducedDim(sce,"tsne",color_by="hto_classification_global")

## ----scds, eval = TRUE, echo = TRUE, message = FALSE, warning = FALSE---------
#- Annotate doublet using co-expression based doublet scoring:
sce = cxds(sce,retRes = TRUE)
sce = bcds(sce,retRes = TRUE,verb=TRUE)
sce = cxds_bcds_hybrid(sce)
par(mfcol=c(1,3))
boxplot(sce$cxds_score   ~ sce$doublet_true_labels, main="cxds")
boxplot(sce$bcds_score   ~ sce$doublet_true_labels, main="bcds")
boxplot(sce$hybrid_score ~ sce$doublet_true_labels, main="hybrid")


## ----pairplot, eval = FALSE, echo = TRUE, message = FALSE, warning = FALSE----
#  scds =
#  top3 = metadata(sce)$cxds$topPairs[1:3,]
#  rs   = rownames(sce)
#  hb   = rowData(sce)$cxds_hvg_bool
#  ho   = rowData(sce)$cxds_hvg_ordr[hb]
#  hgs  = rs[ho]
#  
#  l1 =  ggdraw() + draw_text("Pair 1", x = 0.5, y = 0.5)
#  p1 = plotReducedDim(sce,"tsne",color_by=hgs[top3[1,1]])
#  p2 = plotReducedDim(sce,"tsne",color_by=hgs[top3[1,2]])
#  
#  l2 =  ggdraw() + draw_text("Pair 2", x = 0.5, y = 0.5)
#  p3 = plotReducedDim(sce,"tsne",color_by=hgs[top3[2,1]])
#  p4 = plotReducedDim(sce,"tsne",color_by=hgs[top3[2,2]])
#  
#  l3 = ggdraw() + draw_text("Pair 3", x = 0.5, y = 0.5)
#  p5 = plotReducedDim(sce,"tsne",color_by=hgs[top3[3,1]])
#  p6 = plotReducedDim(sce,"tsne",color_by=hgs[top3[3,2]])
#  
#  plot_grid(l1,p1,p2,l2,p3,p4,l3,p5,p6,ncol=3, rel_widths = c(1,2,2))

## ----sessionInfo--------------------------------------------------------------
sessionInfo()


