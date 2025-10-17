## ----load-libs, message = FALSE,  warning = FALSE-----------------------------
library(ggplot2)
library(SPOTlight)
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)

## ----load-sp, message=FALSE---------------------------------------------------
library(TENxVisiumData)
spe <- MouseKidneyCoronal()
# Use symbols instead of Ensembl IDs as feature names
rownames(spe) <- rowData(spe)$symbol

## ----load-sc, message=FALSE---------------------------------------------------
library(TabulaMurisSenisData)
sce <- TabulaMurisSenisDroplet(tissues = "Kidney")$Kidney

## ----explo--------------------------------------------------------------------
table(sce$free_annotation, sce$age)

## ----sub-18m------------------------------------------------------------------
# Keep cells from 18m mice
sce <- sce[, sce$age == "18m"]
# Keep cells with clear cell type annotations
sce <- sce[, !sce$free_annotation %in% c("nan", "CD45")]

## ----lognorm------------------------------------------------------------------
sce <- logNormCounts(sce)

## ----variance-----------------------------------------------------------------
# Get vector indicating which genes are neither ribosomal or mitochondrial
genes <- !grepl(pattern = "^Rp[l|s]|Mt", x = rownames(sce))

dec <- modelGeneVar(sce, subset.row = genes)
plot(dec$mean, dec$total, xlab = "Mean log-expression", ylab = "Variance")
curve(metadata(dec)$trend(x), col = "blue", add = TRUE)

# Get the top 3000 genes.
hvg <- getTopHVGs(dec, n = 3000)

## ----mgs----------------------------------------------------------------------
colLabels(sce) <- colData(sce)$free_annotation

# Compute marker genes
mgs <- scoreMarkers(sce, subset.row = genes)

## ----mgs-df-------------------------------------------------------------------
mgs_fil <- lapply(names(mgs), function(i) {
    x <- mgs[[i]]
    # Filter and keep relevant marker genes, those with AUC > 0.8
    x <- x[x$mean.AUC > 0.8, ]
    # Sort the genes from highest to lowest weight
    x <- x[order(x$mean.AUC, decreasing = TRUE), ]
    # Add gene and cluster id to the dataframe
    x$gene <- rownames(x)
    x$cluster <- i
    data.frame(x)
})
mgs_df <- do.call(rbind, mgs_fil)

## ----downsample---------------------------------------------------------------
# split cell indices by identity
idx <- split(seq(ncol(sce)), sce$free_annotation)
# downsample to at most 20 per identity & subset
# We are using 5 here to speed up the process but set to 75-100 for your real
# life analysis
n_cells <- 5
cs_keep <- lapply(idx, function(i) {
    n <- length(i)
    if (n < n_cells)
        n_cells <- n
    sample(i, n_cells)
})
sce <- sce[, unlist(cs_keep)]

## ----SPOTlight----------------------------------------------------------------
res <- SPOTlight(
    x = sce,
    y = spe,
    groups = as.character(sce$free_annotation),
    mgs = mgs_df,
    hvg = hvg,
    weight_id = "mean.AUC",
    group_id = "cluster",
    gene_id = "gene")

## ----SPOTligh2, eval=FALSE----------------------------------------------------
#  mod_ls <- trainNMF(
#      x = sce,
#      y = spe,
#      groups = sce$type,
#      mgs = mgs,
#      weight_id = "weight",
#      group_id = "type",
#      gene_id = "gene")
#  
#   # Run deconvolution
#  res <- runDeconvolution(
#      x = spe,
#      mod = mod_ls[["mod"]],
#      ref = mod_ls[["topic"]])

## -----------------------------------------------------------------------------
# Extract deconvolution matrix
head(mat <- res$mat)[, seq_len(3)]
# Extract NMF model fit
mod <- res$NMF

## ----plotTopicProfiles1, fig.width=6, fig.height=7----------------------------
plotTopicProfiles(
    x = mod,
    y = sce$free_annotation,
    facet = FALSE,
    min_prop = 0.01,
    ncol = 1) +
    theme(aspect.ratio = 1)

## ----plotTopicProfiles2, fig.width=9, fig.height=6----------------------------
plotTopicProfiles(
    x = mod,
    y = sce$free_annotation,
    facet = TRUE,
    min_prop = 0.01,
    ncol = 6)

## ----basis-dt, message=FALSE, warning=FALSE-----------------------------------
library(NMF)
sign <- basis(mod)
colnames(sign) <- paste0("Topic", seq_len(ncol(sign)))
head(sign)
# This can be dynamically visualized with DT as shown below
# DT::datatable(sign, fillContainer = TRUE, filter = "top")

## ----plotCorrelationMatrix, fig.width=9, fig.height=9-------------------------
plotCorrelationMatrix(mat)

## ----plotInteractions, fig.width=9, fig.height=9------------------------------
plotInteractions(mat, which = "heatmap", metric = "prop")
plotInteractions(mat, which = "heatmap", metric = "jaccard")
plotInteractions(mat, which = "network")

## ----Scatterpie, fig.width=9, fig.height=6------------------------------------
ct <- colnames(mat)
mat[mat < 0.1] <- 0

# Define color palette
# (here we use 'paletteMartin' from the 'colorBlindness' package)
paletteMartin <- c(
    "#000000", "#004949", "#009292", "#ff6db6", "#ffb6db", 
    "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", 
    "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d")

pal <- colorRampPalette(paletteMartin)(length(ct))
names(pal) <- ct

plotSpatialScatterpie(
    x = spe,
    y = mat,
    cell_types = colnames(mat),
    img = FALSE,
    scatterpie_alpha = 1,
    pie_scale = 0.4) +
    scale_fill_manual(
        values = pal,
        breaks = names(pal))

## -----------------------------------------------------------------------------
plotSpatialScatterpie(
    x = spe,
    y = mat,
    cell_types = colnames(mat),
    img = FALSE,
    scatterpie_alpha = 1,
    pie_scale = 0.4, 
    # Rotate the image 90 degrees counterclockwise
    degrees = -90,
    # Pivot the image on its x axis
    axis = "h") +
    scale_fill_manual(
        values = pal,
        breaks = names(pal))


## ----message=FALSE------------------------------------------------------------
spe$res_ss <- res[[2]][colnames(spe)]
xy <- spatialCoords(spe)
spe$x <- xy[, 1]
spe$y <- xy[, 2]
ggcells(spe, aes(x, y, color = res_ss)) +
    geom_point() +
    scale_color_viridis_c() +
    coord_fixed() +
    theme_bw()

## ----session-info-------------------------------------------------------------
sessionInfo()

