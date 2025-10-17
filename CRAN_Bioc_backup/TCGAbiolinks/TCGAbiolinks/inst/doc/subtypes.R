## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(dpi = 300)
knitr::opts_chunk$set(cache=FALSE)

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)

## ----message=FALSE, warning=FALSE---------------------------------------------
subtypes <- PanCancerAtlas_subtypes()
DT::datatable(
    data = subtypes,
    filter = 'top',
    options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
    rownames = FALSE
)

## ----eval = TRUE--------------------------------------------------------------
lgg.gbm.subtype <- TCGAquery_subtype(tumor = "lgg")

## ----eval = TRUE, echo = FALSE------------------------------------------------
datatable(
    lgg.gbm.subtype[1:10,],
    caption = "Table with LGG molecular subtypes from TCGAquery_subtype",
    filter = 'top',
    options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
    rownames = FALSE
)

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

