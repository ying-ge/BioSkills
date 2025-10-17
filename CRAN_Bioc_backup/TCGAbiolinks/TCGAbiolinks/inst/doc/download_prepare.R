## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(progress = FALSE)

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)

## ----results = 'hide', message=FALSE, warning=FALSE, eval=FALSE---------------
#  # Gene expression aligned against hg38
#  query <- GDCquery(
#      project = "TCGA-GBM",
#      data.category = "Transcriptome Profiling",
#      data.type = "Gene Expression Quantification",
#      workflow.type = "STAR - Counts",
#      barcode = c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01")
#  )
#  GDCdownload(query = query)
#  data <- GDCprepare(query = query)

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
data <- gbm.exp.harmonized

## ----message=FALSE, warning=FALSE---------------------------------------------
datatable(
    as.data.frame(colData(data)), 
    options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
    rownames = FALSE
)

datatable(
    assay(data)[1:20,], 
    options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
    rownames = TRUE
)

rowRanges(data)

## ----eval = FALSE-------------------------------------------------------------
#  query <- GDCquery(
#      project = "TCGA-ACC",
#      data.category = "Copy Number Variation",
#      data.type = "Copy Number Segment",
#      barcode = c( "TCGA-OR-A5KU-01A-11D-A29H-01", "TCGA-OR-A5JK-01A-11D-A29H-01")
#  )
#  GDCdownload(query)
#  data <- GDCprepare(query)

## ----eval = FALSE-------------------------------------------------------------
#  query <- GDCquery(
#      project = "TCGA-ACC",
#      data.category = "Copy Number Variation",
#      data.type = "Gene Level Copy Number",
#      access = "open"
#  )
#  GDCdownload(query)
#  data <- GDCprepare(query)

## ----eval = FALSE-------------------------------------------------------------
#  query <- GDCquery(
#      project = "TCGA-ACC",
#      data.category = "Copy Number Variation",
#      data.type = "Allele-specific Copy Number Segment",
#      access = "open"
#  )
#  GDCdownload(query)
#  data <- GDCprepare(query)

## ----eval = FALSE-------------------------------------------------------------
#  query <- GDCquery(
#      project = "TCGA-ACC",
#      data.category = "Copy Number Variation",
#      data.type = "Masked Copy Number Segment",
#      access = "open"
#  )
#  GDCdownload(query)
#  data <- GDCprepare(query)

## ----eval = FALSE-------------------------------------------------------------
#  # mRNA pipeline: https://gdc-docs.nci.nih.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
#  query.exp.hg38 <- GDCquery(
#      project = "TCGA-GBM",
#      data.category = "Transcriptome Profiling",
#      data.type = "Gene Expression Quantification",
#      workflow.type = "STAR - Counts",
#      barcode =  c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01")
#  )
#  GDCdownload(query.exp.hg38)
#  expdat <- GDCprepare(
#      query = query.exp.hg38,
#      save = TRUE,
#      save.filename = "exp.rda"
#  )

## ----eval = FALSE-------------------------------------------------------------
#  library(TCGAbiolinks)
#  query.mirna <- GDCquery(
#      project = "TARGET-AML",
#      experimental.strategy = "miRNA-Seq",
#      data.category = "Transcriptome Profiling",
#      barcode = c("TARGET-20-PATDNN","TARGET-20-PAPUNR"),
#      data.type = "miRNA Expression Quantification"
#  )
#  GDCdownload(query.mirna)
#  mirna <- GDCprepare(
#      query = query.mirna,
#      save = TRUE,
#      save.filename = "mirna.rda"
#  )
#  

## ----eval = FALSE-------------------------------------------------------------
#  
#  query.isoform <- GDCquery(
#      project = "TARGET-AML",
#      experimental.strategy = "miRNA-Seq",
#      data.category = "Transcriptome Profiling",
#      barcode = c("TARGET-20-PATDNN","TARGET-20-PAPUNR"),
#      data.type = "Isoform Expression Quantification"
#  )
#  GDCdownload(query.isoform)
#  
#  isoform <- GDCprepare(
#      query = query.isoform,
#      save = TRUE,
#      save.filename = "mirna-isoform.rda"
#  )

## ----eval = FALSE-------------------------------------------------------------
#  query_met.hg38 <- GDCquery(
#      project = "TCGA-BRCA",
#      data.category = "DNA Methylation",
#      data.type = "Methylation Beta Value",
#      platform = "Illumina Human Methylation 27",
#      barcode = c("TCGA-B6-A0IM")
#  )
#  GDCdownload(query_met.hg38)
#  data.hg38 <- GDCprepare(query_met.hg38)
#  
#  query_met.hg38 <- GDCquery(
#      project= "TCGA-LGG",
#      data.category = "DNA Methylation",
#      data.type = "Methylation Beta Value",
#      platform = "Illumina Human Methylation 450",
#      barcode = c("TCGA-HT-8111-01A-11D-2399-05","TCGA-HT-A5R5-01A-11D-A28N-05")
#  )
#  GDCdownload(query_met.hg38)
#  data.hg38 <- GDCprepare(query_met.hg38)
#  
#  
#  query_met.hg38 <- GDCquery(
#      project= "HCMI-CMDC",
#      data.category = "DNA Methylation",
#      data.type = "Methylation Beta Value",
#      platform = "Illumina Methylation Epic",
#      barcode = c("HCM-BROD-0045")
#  )
#  GDCdownload(query_met.hg38)
#  data.hg38 <- GDCprepare(query_met.hg38)

## ----eval = FALSE-------------------------------------------------------------
#  query <- GDCquery(
#      project = "TCGA-BRCA",
#      data.category = "DNA Methylation",
#      data.type = "Masked Intensities",
#      platform = "Illumina Human Methylation 27"
#  )
#  GDCdownload(query, files.per.chunk=10)
#  betas <- GDCprepare(query)
#  
#  query <- GDCquery(
#      project = "HCMI-CMDC",
#      data.category = "DNA Methylation",
#      data.type = "Masked Intensities",
#      platform = "Illumina Methylation Epic"
#  )
#  GDCdownload(query, files.per.chunk = 10)
#  betas <- GDCprepare(query)
#  
#  
#  query <- GDCquery(
#      project = "CPTAC-3",
#      data.category = "DNA Methylation",
#      data.type = "Masked Intensities",
#      platform = "Illumina Methylation Epic"
#  )
#  GDCdownload(query, files.per.chunk=10)
#  betas <- GDCprepare(query)
#  
#  query <- GDCquery(
#      project = "TCGA-BRCA",
#      data.category = "DNA Methylation",
#      data.type = "Masked Intensities",
#      platform = "Illumina Methylation Epic"
#  )
#  GDCdownload(query, files.per.chunk = 10)
#  betas <- GDCprepare(query)
#  
#  

## ----eval = FALSE-------------------------------------------------------------
#  query.rppa <- GDCquery(
#      project = "TCGA-ESCA",
#      data.category = "Proteome Profiling",
#      data.type = "Protein Expression Quantification"
#  )
#  GDCdownload(query.rppa)
#  rppa <- GDCprepare(query.rppa)

## ----eval = FALSE-------------------------------------------------------------
#  
#  query <- GDCquery(
#      project = "TCGA-COAD",
#      data.category = "Clinical",
#      data.type = "Clinical Supplement",
#      data.format = "BCR XML",
#      barcode = "TCGA-A6-5664"
#  )
#  GDCdownload(query)
#  drug <- GDCprepare_clinic(query,"drug")
#  
#  query <- GDCquery(
#      project = "TCGA-COAD",
#      data.category = "Clinical",
#      data.type = "Clinical Supplement",
#      data.format = "BCR OMF XML",
#      barcode = "TCGA-AD-6964"
#  )
#  GDCdownload(query)
#  
#  
#  query <- GDCquery(
#      project = "TCGA-ACC",
#      data.category = "Clinical",
#      data.type = "Clinical Supplement",
#      data.format = "BCR Biotab"
#  )
#  GDCdownload(query)
#  clinical.BCRtab.all <- GDCprepare(query)
#  names(clinical.BCRtab.all)
#  
#  query <- GDCquery(
#      project = "TCGA-ACC",
#      data.category = "Clinical",
#      data.type = "Clinical Supplement",
#      data.format = "BCR Biotab",
#      file.type = "radiation"
#  )
#  GDCdownload(query)
#  clinical.BCRtab.radiation <- GDCprepare(query)

## ----eval = FALSE-------------------------------------------------------------
#  query <- GDCquery(
#      project = "TCGA-HNSC",
#      data.category = "Simple Nucleotide Variation",
#      data.type = "Masked Somatic Mutation",
#      access = "open"
#  )
#  GDCdownload(query)
#  maf <- GDCprepare(query)

## ----eval = FALSE-------------------------------------------------------------
#  query.sc.analysis <- GDCquery(
#      project = "CPTAC-3",
#      data.category = "Transcriptome Profiling",
#      access = "open",
#      data.type = "Single Cell Analysis",
#      data.format =  "TSV"
#  )
#  GDCdownload(query.sc.analysis)
#  Single.Cell.Analysis.list <- GDCprepare(query.sc.analysis)

## ----eval = FALSE,include=F---------------------------------------------------
#  query.hdF5 <- GDCquery(
#      project = "CPTAC-3",
#      data.category = "Transcriptome Profiling",
#      access = "open",
#      data.type = "Single Cell Analysis",
#      barcode = c("CPT0167860015","CPT0206880004"),
#      data.format =  "HDF5"
#  )
#  GDCdownload(query.hdF5)
#  df.HDF5 <- GDCprepare(query.hdF5)

## ----eval = FALSE-------------------------------------------------------------
#  query.raw.counts <- GDCquery(
#      project = "CPTAC-3",
#      data.category = "Transcriptome Profiling",
#      access = "open",
#      data.type = "Gene Expression Quantification",
#      barcode = c("CPT0167860015","CPT0206880004"),
#      workflow.type = "CellRanger - 10x Raw Counts"
#  )
#  GDCdownload(query.raw.counts)
#  raw.counts.list <- GDCprepare(query.raw.counts)

## ----eval = FALSE-------------------------------------------------------------
#  query.filtered.counts <- GDCquery(
#      project = "CPTAC-3",
#      data.category = "Transcriptome Profiling",
#      access = "open",
#      data.type = "Gene Expression Quantification",
#      barcode = c("CPT0167860015","CPT0206880004"),
#      workflow.type = "CellRanger - 10x Filtered Counts"
#  )
#  GDCdownload(query.filtered.counts)
#  filtered.counts.list <- GDCprepare(query.filtered.counts)

## ----eval = FALSE-------------------------------------------------------------
#  query.sc.dea <- GDCquery(
#      project = "CPTAC-3",
#      data.category = "Transcriptome Profiling",
#      access = "open",
#      data.type = "Differential Gene Expression",
#      barcode = c("CPT0167860015","CPT0206880004"),
#      workflow.type = "Seurat - 10x Chromium"
#  )
#  GDCdownload(query.sc.dea)
#  sc.dea.list <- GDCprepare(query.sc.dea)

