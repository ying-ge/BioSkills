## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(dpi = 300)
knitr::opts_chunk$set(cache = FALSE)

## ----echo = FALSE,hide=TRUE, message=FALSE,warning=FALSE----------------------
library(TCGAbiolinks)

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
library(SummarizedExperiment)
library(dplyr)
library(DT)

## ----eval = TRUE, message = FALSE, results = "hide"---------------------------
# Selecting TCGA breast cancer (10 samples) for example stored in dataBRCA
dataNorm <- TCGAanalyze_Normalization(
    tabDF = dataBRCA, 
    geneInfo =  geneInfo
)

# quantile filter of genes
dataFilt <- TCGAanalyze_Filtering(
  tabDF = dataNorm,
  method = "quantile",
  qnt.cut =  0.25
)

data(SC_PCBC_stemSig)
Stemness_score <- TCGAanalyze_Stemness(
  stemSig = SC_PCBC_stemSig,
  dataGE = dataFilt
)
data(ECTO_PCBC_stemSig)
ECTO_score <- TCGAanalyze_Stemness(
  stemSig = ECTO_PCBC_stemSig,
  dataGE = dataFilt,
  colname.score = "ECTO_PCBC_stem_score"
)

data(MESO_PCBC_stemSig)
MESO_score <- TCGAanalyze_Stemness(
  stemSig = MESO_PCBC_stemSig,
  dataGE = dataFilt,
  colname.score = "MESO_PCBC_stem_score"
)

## ----eval = T-----------------------------------------------------------------
head(Stemness_score)
head(ECTO_score)
head(MESO_score)

