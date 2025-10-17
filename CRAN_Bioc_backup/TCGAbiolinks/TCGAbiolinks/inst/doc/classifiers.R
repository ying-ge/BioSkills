## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(dpi = 300)
knitr::opts_chunk$set(cache = FALSE)

## ----echo = FALSE,hide=TRUE, message=FALSE,warning=FALSE----------------------
library(TCGAbiolinks)

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
library(SummarizedExperiment)
library(dplyr)
library(DT)

## ----eval = FALSE, message = FALSE, results = "hide"--------------------------
#  query <- GDCquery(
#      project = "TCGA-GBM",
#      data.category = "DNA Methylation",
#      barcode = c("TCGA-06-0122","TCGA-14-1456"),
#      platform = "Illumina Human Methylation 27",
#      data.type = "Methylation Beta Value"
#  )
#  GDCdownload(query)
#  dnam <- GDCprepare(query)

## ----eval = FALSE-------------------------------------------------------------
#  assay(dnam)[1:5,1:2]

## ----eval = FALSE-------------------------------------------------------------
#  classification <- gliomaClassifier(dnam)

## ----eval = FALSE-------------------------------------------------------------
#  names(classification)
#  classification$final.classification
#  classification$model.classifications
#  classification$model.probabilities

## -----------------------------------------------------------------------------
TCGAquery_subtype("GBM") %>%
 dplyr::filter(patient %in% c("TCGA-06-0122","TCGA-14-1456")) %>%
 dplyr::select("patient","Supervised.DNA.Methylation.Cluster")

