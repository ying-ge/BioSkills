## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(progress = FALSE)

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)

## ----results = 'hide', echo=TRUE, message=FALSE, warning=FALSE,eval=F---------
#  query <- GDCquery(
#      project = "TCGA-CHOL",
#      data.category = "Simple Nucleotide Variation",
#      access = "open",
#      data.type = "Masked Somatic Mutation",
#      workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
#  )
#  GDCdownload(query)
#  maf <- GDCprepare(query)

## ----results = 'hide', echo=TRUE, message=FALSE, warning=FALSE,eval=T,include=F----
maf <- chol_maf@data

## ----echo = TRUE, message = FALSE, warning = FALSE----------------------------
# Only first 50 to make render faster
datatable(maf[1:20,],
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

## ----results = 'hide', echo=TRUE, message=FALSE, warning=FALSE,eval=FALSE-----
#  maf <- getMC3MAF()

## ----results = "hide",echo = TRUE, message = FALSE, warning = FALSE, eval=FALSE----
#  library(maftools)
#  library(dplyr)
#  query <- GDCquery(
#      project = "TCGA-CHOL",
#      data.category = "Simple Nucleotide Variation",
#      access = "open",
#      data.type = "Masked Somatic Mutation",
#      workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
#  )
#  GDCdownload(query)
#  maf <- GDCprepare(query)
#  
#  maf <- maf %>% maftools::read.maf

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
library(maftools)
library(dplyr)
maf <- chol_maf

## ----results = "hide",echo = TRUE, message = FALSE, warning = FALSE-----------
datatable(getSampleSummary(maf),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)

## ----echo = TRUE, message = FALSE,eval = FALSE, warning = FALSE---------------
#  oncoplot(maf = maf, top = 10, removeNonMutated = TRUE)
#  titv = titv(maf = maf, plot = FALSE, useSyn = TRUE)
#  #plot titv summary
#  plotTiTv(res = titv)

