## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(cache = TRUE)

## ----warning = FALSE, message = FALSE-----------------------------------------
# load required packages
library(CATALYST)
library(cowplot)
library(flowCore)
library(diffcyt)
library(scater)
library(SingleCellExperiment)

## ----load-data----------------------------------------------------------------
# load example data
data(PBMC_fs, PBMC_panel, PBMC_md)
PBMC_fs
head(PBMC_panel)
head(PBMC_md)

## ----eval=FALSE---------------------------------------------------------------
#  # download exemplary set of FCS files
#  url <- "http://imlspenticton.uzh.ch/robinson_lab/cytofWorkflow"
#  zip <- "PBMC8_fcs_files.zip"
#  download.file(paste0(url, "/", zip), destfile = zip, mode = "wb")
#  unzip(zip)
#  
#  # read in FCS files as flowSet
#  fcs <- list.files(pattern = ".fcs$")
#  fs <- read.flowSet(fcs, transformation = FALSE, truncate_max_range = FALSE)

## ----prepData-----------------------------------------------------------------
(sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md))

## ----eval=FALSE---------------------------------------------------------------
#  # alter panel column names
#  panel2 <- PBMC_panel
#  colnames(panel2)[1:2] <- c("channel_name", "marker")
#  
#  # alter metadata column names & add 2nd condition
#  md2 <- PBMC_md
#  colnames(md2) <- c("file", "sampleID", "cond1", "patientID")
#  md2$cond2 <- rep(c("A", "B"), 4)
#  
#  # construct SCE
#  prepData(PBMC_fs, panel2, md2,
#      panel_cols = list(channel = "channel_name", antigen = "marker"),
#      md_cols = list(file = "file", id = "sampleID",
#          factors = c("cond1", "cond2", "patientID")))

## ----cluster------------------------------------------------------------------
sce <- cluster(sce, features = "type", 
    xdim = 10, ydim = 10, maxK = 20, 
    verbose = FALSE, seed = 1)       

## ----mergeClusters------------------------------------------------------------
data(merging_table)
head(merging_table)
sce <- mergeClusters(sce, k = "meta20", table = merging_table, id = "merging1")
head(cluster_codes(sce))[, seq_len(10)]

## ----delta-area, fig.width = 5, fig.height = 3--------------------------------
# access & render delta area plot
# (equivalent to metadata(sce)$delta_area)
delta_area(sce)

## ----plotCounts-1, fig.width = 5, fig.height = 3------------------------------
plotCounts(sce, 
    group_by = "sample_id", 
    color_by = "condition")

## ----plotCounts-2, fig.width = 4, fig.height = 3------------------------------
plotCounts(sce, 
    prop = TRUE,
    group_by = "condition", 
    color_by = "patient_id")

## ----pbMDS-1, fig.width = 5---------------------------------------------------
pbMDS(sce, shape_by = "patient_id", size_by = FALSE)

## ----pbMDS-2, fig.width = 5---------------------------------------------------
pbMDS(sce, by = "both", k = "meta12", 
    shape_by = "condition", size_by = TRUE)

## ----clrDR-1, fig.width = 5---------------------------------------------------
clrDR(sce, by = "sample_id", k = "meta12")

## ----clrDR-2, fig.width = 6---------------------------------------------------
clrDR(sce, by = "cluster_id", arrow_col = "condition", size_by = FALSE)

## ----plotExprHeatmap-sample, fig.width = 6.5, fig.height = 3.5----------------
# scale each marker between 0 and 1 
# (after aggregation & without trimming)
plotExprHeatmap(sce, features = "state", 
    scale = "last", q = 0, bars = FALSE)

## ----plotExprHeatmap-cluster, fig.width = 6, fig.height = 3.5-----------------
# medians of scaled & trimmed type-marker expressions by cluster
plotExprHeatmap(sce, features = "type",
    by = "cluster_id", k = "meta12", m = "meta8",
    scale = "first", q = 0.01, perc = TRUE, col_dend = FALSE)

## ----plotExprHeatmap-both, fig.width = 4, fig.height = 3.5--------------------
# raw (not scaled, not trimmed) 
# median expression by cluster-sample
plotExprHeatmap(sce, features = "pS6", by = "both", k = "meta8", 
    scale = "never", col_clust = FALSE, row_anno = FALSE, bars = FALSE)

## ----plotMedExprs-1, fig.width = 12, fig.height = 5---------------------------
plotPbExprs(sce, k = "meta8", facet_by = "cluster_id", ncol = 4)

## ----plotMedExprs-2, fig.width = 12, fig.height = 4---------------------------
plotPbExprs(sce, facet_by = "antigen", ncol = 7)

## ----plotMedExprs-3, fig.width = 12, fig.height = 4---------------------------
plotPbExprs(sce, k = "meta10", features = "type", 
  group_by = "cluster_id", color_by = "sample_id", 
  size_by = TRUE, geom = "points", jitter = FALSE, ncol = 5)

## ----plotMedExprs-4, fig.width = 12, fig.height = 4---------------------------
plotPbExprs(sce, k = "meta6", features = "state",  
    group_by = "cluster_id", color_by = "condition", ncol = 7)

## ----plotClusterExprs, message = FALSE, fig.width = 12, fig.height = 8--------
plotClusterExprs(sce, k = "meta8", features = "type")

## ----plotAbundances-1, fig.width = 5, fig.height = 3.5------------------------
plotAbundances(sce, k = "meta12", by = "sample_id", group_by = "condition")

## ----plotAbundances-2, fig.width = 6, fig.height = 4--------------------------
plotAbundances(sce, k = "meta8", by = "cluster_id", 
    group_by = "condition", shape_by = "patient_id")

## ----plotFreqHeatmap-complete, fig.width = 6.5, fig.height = 4.5--------------
# complete example
plotFreqHeatmap(sce, 
    k = "meta8", m = "meta5",
    hm_pal = rev(hcl.colors(10, "RdBu")),
    k_pal = hcl.colors(7, "Zissou 1"),
    m_pal = hcl.colors(4, "Temps"),
    bars = TRUE, perc = TRUE)

## ----plotFreqHeatmap-minimal, fig.width = 3.5, fig.height = 4-----------------
# minimal example
plotFreqHeatmap(sce, k = "meta10", 
    normalize = FALSE, bars = FALSE,
    row_anno = FALSE, col_anno = FALSE,
    row_clust = FALSE, col_clust = FALSE,
    hm_pal = c("grey95", "black"))

## ----plotMultiHeatmap-1, fig.width = 7.5, fig.height = 4----------------------
# both, median type- & state-marker expressions
plotMultiHeatmap(sce, 
    hm1 = "type", hm2 = "state", 
    k = "meta12", m = "meta8",
    col_dend = c(FALSE, TRUE))

## ----plotMultiHeatmap-2, fig.width = 7.5, fig.height = 3----------------------
# 1st: CDx markers by cluster; 
# 2nd: population frequencies by sample
cdx <- grep("CD", rownames(sce), value = TRUE)
plotMultiHeatmap(sce, k = "meta6",
    hm1 = cdx, hm2 = "abundances", 
    bars = TRUE, perc = TRUE, row_anno = FALSE)

## ----plotMultiHeatmap-3, fig.width = 9, fig.height = 4------------------------
# plot selected markers side-by-side;
# omit left-hand side heatmap
plotMultiHeatmap(sce, 
    k = "meta8", scale = "never",
    hm1 = FALSE, hm2 = c("pS6", "pp38", "pBtk"),
    row_anno = FALSE, col_clust = FALSE,
    hm2_pal = c("grey95", "black"))

## ----runDR--------------------------------------------------------------------
set.seed(1601)
sce <- runDR(sce, dr = "UMAP", cells = 500, features = "type")

## ----runUMAP-scater, eval = FALSE---------------------------------------------
#  sce <- runUMAP(sce, exprs_values = "exprs")

## -----------------------------------------------------------------------------
# view & access DRs
reducedDimNames(sce)
head(reducedDim(sce, "UMAP"))

## ----plotDR-1, fig.height = 5-------------------------------------------------
# color by marker expression & split by condition
plotDR(sce, color_by = c("pS6", "pNFkB"), facet_by = "condition")

## ----plotDR-2, fig.height = 5-------------------------------------------------
# color by 8 metaclusters & split by sample ID
p <- plotDR(sce, color_by = "meta8", facet_by = "sample_id")
p$facet$params$ncol <- 4; p

## ----filterSCE, fig.height = 3------------------------------------------------
u <- filterSCE(sce, patient_id == "Patient1")
table(u$sample_id)

u <- filterSCE(sce, k = "meta8",
    cluster_id %in% c(1, 3, 8))
plot_grid(
    plotDR(sce, color_by = "meta8"),
    plotDR(u, color_by = "meta8"))

## ----diffcyt, message = FALSE, warning = FALSE, fig.show = "hide"-------------
# create design & contrast matrix
design <- createDesignMatrix(ei(sce), cols_design = "condition")
contrast <- createContrast(c(0, 1))

# test for
# - differential abundance (DA) of clusters
# - differential states (DS) within clusters
res_DA <- diffcyt(sce, clustering_to_use = "meta10",
    analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
    design = design, contrast = contrast, verbose = FALSE)
res_DS <- diffcyt(sce, clustering_to_use = "meta10",
    analysis_type = "DS", method_DS = "diffcyt-DS-limma",
    design = design, contrast = contrast, verbose = FALSE)

# extract result tables
tbl_DA <- rowData(res_DA$res)
tbl_DS <- rowData(res_DS$res)

## ----plotDiffHeatmap-da, fig.width = 5, fig.height = 4------------------------
plotDiffHeatmap(sce, tbl_DA, all = TRUE, fdr = 0.05)

## ----plotDiffHeatmap-ds, fig.width = 5, fig.height = 4------------------------
plotDiffHeatmap(sce, tbl_DS, fdr = 0.05, 
    sort_by = "lfc", col_anno = "condition")

## ----plotDiffHeatmap-filter-1, fig.width = 4.5, fig.height = 3.25-------------
# include all results for selected marker
plotDiffHeatmap(sce["pp38", ], tbl_DS, all = TRUE, col_anno = FALSE)

## ----plotDiffHeatmap-filter-2, fig.width = 5, fig.height = 4.25---------------
# include all results for selected cluster
k <- metadata(res_DS$res)$clustering_name
sub <- filterSCE(sce, cluster_id == 8, k = k)
plotDiffHeatmap(sub, tbl_DS, all = TRUE, normalize = FALSE)

## ----plotDiffHeatmap-pals, fig.width = 4, fig.height = 3.5--------------------
plotDiffHeatmap(sce, tbl_DA, all = TRUE, col_anno = FALSE,
    hm_pal = c("gold", "white", "navy"),
    fdr_pal = c("grey90", "grey50"),
    lfc_pal = c("red3", "grey90", "green3"))

## ----sce2fcs-split, message = FALSE-------------------------------------------
# store final clustering in cell metadata
sce$mm <- cluster_ids(sce, "merging1")
# convert to 'flowSet' with one frame per cluster 
(fs <- sce2fcs(sce, split_by = "mm"))
# split check: number of cells per barcode ID
# equals number of cells in each 'flowFrame'
all(c(fsApply(fs, nrow)) == table(sce$mm))
# store identifiers (= cluster names)
(ids <- c(fsApply(fs, identifier)))

## ----sce2fcs-write, eval = FALSE----------------------------------------------
#  for (id in ids) {
#      # subset 'flowFrame' for cluster 'id'
#      ff <- fs[[id]]
#      # specify output name that includes ID
#      fn <- sprintf("manuel_merging_%s.fcs", id)
#      # construct output path
#      fn <- file.path("...", fn)
#      # write frame to FCS
#      write.FCS(ff, fn)
#  }

## ----other-clusterings, message = FALSE, warning = FALSE----------------------
# subset type-marker expression matrix
es <- assay(sce, "exprs")
es <- es[type_markers(sce), ]

# run clustering method X
# (here, we just split the cells into 
# equal chunks according to CD33 expression)
cs <- split(seq_len(ncol(sce)), cut(es["CD33", ], nk <- 10))
kids <- lapply(seq_len(nk), function(i) {
    rep(i, length(cs[[i]]))
})
kids <- factor(unlist(kids))

# store cluster IDs in cell metadata & codes in metadata
foo <- sce
foo$cluster_id[unlist(cs)] <- unlist(kids)
metadata(foo)$cluster_codes <- data.frame(
    custom = factor(levels(kids), levels = levels(kids)))

# tabulate cluster assignments
table(cluster_ids(foo, "custom"))

## ----plotMedExprs-custom, fig.width = 8, fig.height = 2-----------------------
p <- plotMedExprs(sce, k = "meta4", facet_by = "cluster_id")
# facetting layout is 2x2; plot all side-by-side instead
p$facet$params$nrow <- 1
# remove points
p$layers <- p$layers[-1]
# overwrite default colors
p <- p + scale_color_manual(values = c("royalblue", "orange"))
# remove x-axis title, change angle & decrease size of labels
(p + labs(x = NULL) + theme(axis.text.x = element_text(angle = 90, size = 8)))

## ----plotMultiHeatmap-custom, fig.width = 6, fig.height = 3-------------------
plotMultiHeatmap(sce,
    k = "meta8", 
    m = "meta4",
    hm2 = "abundances",
    # include all dendrograms
    row_dend = TRUE, 
    col_dend = TRUE, 
    # exclude sample annotations
    col_anno = FALSE,
    # primary & merging cluster palettes
    k_pal = hcl.colors(8, "Vik"),     
    m_pal = hcl.colors(4, "Tropic"), 
    # 1st & 2nd heatmap coloring
    hm1_pal = c("grey95", "blue"),  
    hm2_pal = c("grey95", "red3"))

## ----plotExprHeatmap-minimal, fig.width = 8, fig.height = 3-------------------
# minimal heatmap
plotExprHeatmap(sce,
    row_anno = FALSE,   # don't annotate samples
    row_clust = FALSE,  # keep samples in original order
    col_clust = FALSE,  # keep markers in original order
    bin_anno = FALSE,   # don't annotate bins
    bars = FALSE,       # don't include sample sizes
    scale = "last",     # aggregate, then scale
    hm_pal = hcl.colors(10, "YlGnBu", rev = TRUE))

## ----plotExprHeatmap-complete, fig.width = 12, fig.height = 4-----------------
# complete heatmap
plotExprHeatmap(sce, row_anno = TRUE,   # annotate samples
    row_clust = TRUE, col_clust = TRUE, # cluster samples/markers
    row_dend = TRUE, col_dend = TRUE,   # include dendrograms
    bin_anno = TRUE,          # annotate bins with value
    bars = TRUE, perc = TRUE, # include barplot of sample sizes
    hm_pal = c("grey95", "orange"))

## ----combine-heatmaps-1, fig.width = 14, fig.height = 4.8---------------------
# specify clustering to aggregate by
k <- "meta11"

# median type-marker expression by cluster
p1 <- plotExprHeatmap(sce, features = "type",
    by = "cluster_id", k = k, m = "meta7")

# median state-marker expression by cluster
p2 <- plotExprHeatmap(sce, features = "state",
    by = "cluster_id", k = k, row_anno = FALSE)

# relative cluster abundances by sample
p3 <- plotFreqHeatmap(sce, k = k, perc = TRUE,
    row_anno = FALSE, col_clust = FALSE)

# make legend titles unique
p1@name <- p1@matrix_color_mapping@name <- "type"
p2@name <- p2@matrix_color_mapping@name <- "state"

p1 + p2 + p3

## ----combine-heatmaps-2, fig.width = 14, fig.height = 4-----------------------
# specify clustering to aggregate by
k <- "meta9"

# relative cluster abundances by sample
p <- plotFreqHeatmap(sce, k = k, 
    bars = FALSE, hm_pal = c("white", "black"),
    row_anno = FALSE, col_clust = FALSE, col_anno = FALSE)

# specify unique coloring
cs <- c(pp38 = "maroon", pBtk = "green4")

# median expression of selected markers by cluster-sample
for (f in names(cs)) {
    q <- plotExprHeatmap(sce, features = f, 
        by = "both", k = k, scale = "never",
        row_anno = FALSE, col_clust = FALSE, 
        hm_pal = c("white", cs[f]))
    # make identifier & legend title unique
    q@name <- q@matrix_color_mapping@name <- f
    # remove column annotation names
    for (i in seq_along(q@top_annotation@anno_list))
        q@top_annotation@anno_list[[i]]@name_param$show <- FALSE
    # remove redundant sample names
    q@column_names_param$show <- FALSE
    p <- p + q
}

# add heatmap of median expression across all features
p + plotExprHeatmap(sce, features = NULL,
    by = "cluster_id", k = k, row_anno = FALSE)

## ----session-info-------------------------------------------------------------
sessionInfo()

