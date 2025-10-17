## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(progress = FALSE)

## ----message=FALSE, warning=FALSE, eval=FALSE---------------------------------
#  if (!requireNamespace("BiocManager", quietly=TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("TCGAbiolinks")

## ----message=FALSE, warning=FALSE, eval=FALSE---------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
#  BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")

## ----message=FALSE, warning=FALSE, include=TRUE-------------------------------
library(TCGAbiolinks)
library(dplyr)
library(DT)

## -----------------------------------------------------------------------------
version
packageVersion("TCGAbiolinks")

