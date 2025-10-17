## ----style, echo = FALSE, results = 'asis'------------------------------------
BiocStyle::markdown(css.files = c('custom.css'))

## ----echo=FALSE---------------------------------------------------------------
suppressPackageStartupMessages({
  library(tRNA)
  library(Structstrings)
})
data("gr", package = "tRNA")

## ----eval=FALSE---------------------------------------------------------------
#  library(tRNA)
#  library(Structstrings)
#  data("gr", package = "tRNA")

## -----------------------------------------------------------------------------
# just get the coordinates of the anticodonloop
gettRNAstructureGRanges(gr, structure = "anticodonLoop")
gettRNAstructureSeqs(gr, joinFeatures = TRUE, structure = "anticodonLoop")

## -----------------------------------------------------------------------------
seqs <- gettRNAstructureSeqs(gr[1L:10L], joinCompletely = TRUE)
seqs
# getting the tRNA structure boundaries
metadata(seqs)[["tRNA_structures"]]

## ----echo=TRUE, results="hide"------------------------------------------------
gr[hasAcceptorStem(gr, unpaired = TRUE)]
# mismatches and bulged are subsets of unpaired
gr[hasAcceptorStem(gr, mismatches = TRUE)]
gr[hasAcceptorStem(gr, bulged = TRUE)]
# combination of different structure parameters
gr[hasAcceptorStem(gr, mismatches = TRUE) & 
     hasDloop(gr, length = 8L)]

## -----------------------------------------------------------------------------
# load tRNA data for E. coli and H. sapiens
data("gr_eco", package = "tRNA")
data("gr_human", package = "tRNA")

# get summary plots
grl <- GRangesList(Sce = gr,
                   Hsa = gr_human,
                   Eco = gr_eco)
plots <- gettRNAFeaturePlots(grl)

## ----plot1, fig.cap = "tRNA length."------------------------------------------
plots$length

## ----plot2, fig.cap = "tRNAscan-SE scores."-----------------------------------
plots$tRNAscan_score

## ----plot3, fig.cap = "tRNA GC content."--------------------------------------
plots$gc

## ----plot4, fig.cap = "tRNAs with introns."-----------------------------------
plots$tRNAscan_intron

## ----plot5, fig.cap = "Length of the variable loop."--------------------------
plots$variableLoop_length

## ----eval=FALSE---------------------------------------------------------------
#  # score column will be used
#  plots <- gettRNAFeaturePlots(grl, plotScores = TRUE)

## -----------------------------------------------------------------------------
plots <- gettRNAFeaturePlots(grl,
                             scores = list(runif(length(grl[[1L]]),0L,100L),
                                           runif(length(grl[[2L]]),0L,100L),
                                           runif(length(grl[[3L]]),0L,100L)))

## ----plot6, fig.cap = "tRNA length and score correlation."--------------------
plots$length

## ----plot7, fig.cap = "variable loop length and score correlation."-----------
plots$variableLoop_length

## ----plot8, fig.cap = "Customized plot switching out the point and violin plot into a boxplot."----
plots$length$layers <- plots$length$layers[c(-1L,-2L)]
plots$length + ggplot2::geom_boxplot()

## -----------------------------------------------------------------------------
head(plots$length$data)

## -----------------------------------------------------------------------------
options("tRNA_colour_palette")
options("tRNA_colour_yes")
options("tRNA_colour_no")

## -----------------------------------------------------------------------------
head(gettRNABasePairing(gr)[[1L]])
head(getBasePairing(gr[1L]$tRNA_str)[[1L]])

## -----------------------------------------------------------------------------
gettRNALoopIDs(gr)[[1L]]
getLoopIndices(gr[1L]$tRNA_str)

## -----------------------------------------------------------------------------
sessionInfo()

