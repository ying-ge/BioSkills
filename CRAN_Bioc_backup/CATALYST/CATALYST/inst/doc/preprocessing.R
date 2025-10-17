## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(cache = TRUE)

## ----load-libs, warning = FALSE, message = FALSE------------------------------
# load required packages
library(CATALYST)
library(cowplot)
library(flowCore)
library(ggplot2)
library(SingleCellExperiment)

## ----prepData-----------------------------------------------------------------
data("raw_data")
(sce <- prepData(raw_data))
# view number of events per sample
table(sce$sample_id)
# view non-mass channels
names(int_colData(sce))

## ----normCytof, message = FALSE, fig.width = 8, fig.height = 6----------------
# construct SCE
sce <- prepData(raw_data)
# apply normalization; keep raw data
res <- normCytof(sce, beads = "dvs", k = 50, 
  assays = c("counts", "exprs"), overwrite = FALSE)
# check number & percentage of bead / removed events
n <- ncol(sce); ns <- c(ncol(res$beads), ncol(res$removed))
data.frame(
    check.names = FALSE, 
    "#" = c(ns[1], ns[2]), 
    "%" = 100*c(ns[1]/n, ns[2]/n),
    row.names = c("beads", "removed"))
# extract data excluding beads & doublets,
# and including normalized intensitied
sce <- res$data
assayNames(sce)

## ----normCytof-scatter, fig.wide = TRUE, fig.width = 10, fig.height = 2-------
# plot bead vs. dna scatters
res$scatter

## ----normCytof-lines, fig.wide = TRUE, fig.width = 8, fig.height = 4----------
# plot smoothed bead intensities
res$lines

## -----------------------------------------------------------------------------
data(sample_ff)
sample_ff

## -----------------------------------------------------------------------------
data(sample_key)
head(sample_key)

## ----assignPrelim, messages = FALSE-------------------------------------------
sce <- prepData(sample_ff)
(sce <- assignPrelim(sce, sample_key))
# view barcode channels
rownames(sce)[rowData(sce)$is_bc]
# view number of events assigned to each barcode population
table(sce$bc_id)

## ----estCutoffs---------------------------------------------------------------
# estimate separation cutoffs
sce <- estCutoffs(sce)
# view separation cutoff estimates
metadata(sce)$sep_cutoffs

## ----eval = FALSE-------------------------------------------------------------
#  plotYields(sce, which = c(0, "C1"))

## ----plotYields, echo = FALSE, fig.width = 7, fig.height = 3.5----------------
ps <- plotYields(sce, which = c(0, "C1")); ps[[1]]; ps[[2]]

## ----applyCutoffs-------------------------------------------------------------
# use global / population-specific separation cutoff(s)
sce2 <- applyCutoffs(sce)
sce3 <- applyCutoffs(sce, sep_cutoffs = 0.35)

# compare yields before and after applying 
# global / population-specific cutoffs
c(specific = mean(sce2$bc_id != 0),
    global = mean(sce3$bc_id != 0))
# proceed with population-specific filtering
sce <- sce2

## ----eval = FALSE-------------------------------------------------------------
#  # event plots for unassigned events
#  # & barcode population D1
#  plotEvents(sce, which = c(0, "D1"), n = 25)

## ----plotEvents, echo = FALSE, fig.width = 6, fig.height = 3------------------
ps <- plotEvents(sce, which = c(0, "D1"), n = 25); ps[[1]]; ps[[2]]

## ----plotMahal, fig.width = 6, fig.height = 6.5-------------------------------
plotMahal(sce, which = "B3")

## ----computeSpillmat----------------------------------------------------------
# get single-stained control samples
data(ss_exp)

# specify mass channels stained for & debarcode
bc_ms <- c(139, 141:156, 158:176)
sce <- prepData(ss_exp)
sce <- assignPrelim(sce, bc_ms, verbose = FALSE)
sce <- applyCutoffs(estCutoffs(sce))

# compute & extract spillover matrix
sce <- computeSpillmat(sce)
sm <- metadata(sce)$spillover_matrix

# do some sanity checks
chs <- channels(sce)
ss_chs <- chs[rowData(sce)$is_bc]
all(diag(sm[ss_chs, ss_chs]) == 1)
all(sm >= 0 & sm <= 1)

## ----plotSpillmat, fig.width = 6, fig.height = 6------------------------------
plotSpillmat(sce) 

## ----compCytof, message = FALSE, fig.width = 8, fig.height = 3----------------
# construct SCE of multiplexed cells
data(mp_cells)
sce <- prepData(mp_cells)
# compensate using NNLS-method; keep uncompensated data
sce <- compCytof(sce, sm, method = "nnls", overwrite = FALSE)
# visualize data before & after compensation
chs <- c("Er167Di", "Er168Di")
as <- c("exprs", "compexprs")
ps <- lapply(as, function(a) 
    plotScatter(sce, chs, assay = a))
plot_grid(plotlist = ps, nrow = 1)

## ----plotScatter-1, fig.small = TRUE, fig.height = 2.5------------------------
# biscatter of DNA channels colored by cell density
sce <- prepData(raw_data)
chs <- c("DNA1", "DNA2")
plotScatter(sce, chs)

## ----plotScatter-2, fig.width = 6, fig.height = 4-----------------------------
# biscatters for selected CD-channels
sce <- prepData(mp_cells)
chs <- grep("^CD", rownames(sce), value = TRUE)
chs <- sample(chs, 7)
p <- plotScatter(sce, chs)
p$facet$params$ncol <- 3; p

## ----plotScatter-3, message = FALSE, fig.width = 6, fig.height = 3------------
sce <- prepData(sample_ff)
sce <- assignPrelim(sce, sample_key)
# downsample channels & barcode populations
chs <- sample(rownames(sce), 4)
ids <- sample(rownames(sample_key), 3)
sce <- sce[chs, sce$bc_id %in% ids]

# color by factor variable
plotScatter(sce, chs, color_by = "bc_id")

# color by continuous variable
plotScatter(sce, chs, color_by = "delta")

## ----plotScatter-4, message = FALSE, fig.width = 6, fig.height = 4------------
# sample some random group labels
sce$group_id <- sample(c("groupA", "groupB"), ncol(sce), TRUE)

# selected pair of channels; split by barcode & group ID
plotScatter(sce, sample(chs, 2), 
  color_by = "bc_id",
  facet_by = c("bc_id", "group_id"))

## ----plotScatter-5, message = FALSE, fig.width = 6, fig.height = 6------------
# selected CD-channels; split by sample
plotScatter(sce, chs, bins = 50, facet_by = "bc_id")

## ----sce2fcs-split, message = FALSE-------------------------------------------
# run debarcoding
sce <- prepData(sample_ff)
sce <- assignPrelim(sce, sample_key)
sce <- applyCutoffs(estCutoffs(sce))
# exclude unassigned events
sce <- sce[, sce$bc_id != 0]
# convert to 'flowSet' with one frame per sample
(fs <- sce2fcs(sce, split_by = "bc_id"))
# split check: number of cells per barcode ID
# equals number of cells in each 'flowFrame'
all(c(fsApply(fs, nrow)) == table(sce$bc_id))

## ----sce2fcs-write, eval = FALSE----------------------------------------------
#  # get sample identifiers
#  ids <- fsApply(fs, identifier)
#  for (id in ids) {
#      ff <- fs[[id]]                     # subset 'flowFrame'
#      fn <- sprintf("sample_%s.fcs", id) # specify output name that includes ID
#      fn <- file.path("...", fn)         # construct output path
#      write.FCS(ff, fn)                  # write frame to FCS
#  }

## ----sce2fcs-gating, message = FALSE, fig.small = TRUE, fig.height = 2.5------
# load required packages
library(ggcyto)      
library(openCyto)     
library(flowWorkspace) 

# construct 'GatingSet'
sce <- prepData(raw_data) 
ff <- sce2fcs(sce, assay = "exprs")   
gs <- GatingSet(flowSet(ff))

# specify DNA channels
dna_chs <- c("Ir191Di", "Ir193Di")

# apply elliptical gate
gs_add_gating_method(
    gs, alias = "cells", 
    pop = "+", parent = "root",
    dims = paste(dna_chs, collapse = ","),
    gating_method = "flowClust.2d", 
    gating_args = "K=1,q=0.9")

# plot scatter of DNA channels including elliptical gate
ggcyto(gs, 
    aes_string(dna_chs[1], dna_chs[2])) + 
    geom_hex(bins = 128) + 
    geom_gate(data = "cells") +
    facet_null() + ggtitle(NULL) +
    theme(aspect.ratio = 1, 
        panel.grid.minor = element_blank())

## ----session-info-------------------------------------------------------------
sessionInfo()

