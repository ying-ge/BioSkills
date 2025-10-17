## ----setup, echo = FALSE, results = "hide", message = FALSE-------------------
require(knitr)
library(BiocStyle)
opts_chunk$set(error = FALSE, message = FALSE, warning = FALSE)

## ----read---------------------------------------------------------------------
library(zellkonverter)

# Obtaining an example H5AD file.
example_h5ad <- system.file("extdata", "krumsiek11.h5ad",
                            package = "zellkonverter")
readH5AD(example_h5ad)

## ----write--------------------------------------------------------------------
library(scRNAseq)

sce_zeisel <- ZeiselBrainData()
out_path <- tempfile(pattern = ".h5ad")
writeH5AD(sce_zeisel, file = out_path)

## ----convert------------------------------------------------------------------
library(basilisk)
library(scRNAseq)

seger <- SegerstolpePancreasData()
roundtrip <- basiliskRun(fun = function(sce) {
     # Convert SCE to AnnData:
     adata <- SCE2AnnData(sce)

     # Maybe do some work in Python on 'adata':
     # BLAH BLAH BLAH

     # Convert back to an SCE:
     AnnData2SCE(adata)
}, env = zellkonverterAnnDataEnv(), sce = seger)

## ----anndata-deps-------------------------------------------------------------
AnnDataDependencies()

## ----anndata-deps-old---------------------------------------------------------
AnnDataDependencies(version = "0.7.6")

## ----verbose------------------------------------------------------------------
readH5AD(example_h5ad, verbose = TRUE)

## ----verbose-set, eval = FALSE------------------------------------------------
#  # This is not run here
#  setZellkonverterVerbose(TRUE)

## -----------------------------------------------------------------------------
sessionInfo()

