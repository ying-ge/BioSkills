## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(progress = FALSE)

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)

## ----eval = TRUE, echo = FALSE------------------------------------------------
datatable(
    TCGAbiolinks:::getGDCprojects(),
    filter = 'top',
    options = list(scrollX = TRUE, keys = TRUE, pageLength = 10), 
    rownames = FALSE,
    caption = "List of projects"
)

## ----eval = TRUE, echo = FALSE------------------------------------------------
datatable(
    TCGAbiolinks:::getBarcodeDefinition(),
    filter = 'top',
    options = list(scrollX = TRUE, keys = TRUE, pageLength = 10), 
    rownames = FALSE,
    caption = "List sample types"
)

## ----echo=FALSE---------------------------------------------------------------
datatable(
    readr::read_csv("https://docs.google.com/spreadsheets/d/1f98kFdj9mxVDc1dv4xTZdx8iWgUiDYO-qiFJINvmTZs/export?format=csv&gid=2046985454",col_types = readr::cols()),
    filter = 'top',
    options = list(scrollX = TRUE, keys = TRUE, pageLength = 40), 
    rownames = FALSE
)

## ----message=FALSE, warning=FALSE---------------------------------------------
query <- GDCquery(
    project = c("TCGA-GBM", "TCGA-LGG"),
    data.category = "DNA Methylation",
    platform = c("Illumina Human Methylation 450"),
    sample.type = "Recurrent Tumor"
)
datatable(
    getResults(query), 
    filter = 'top',
    options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
    rownames = FALSE
)

## ----message=FALSE, warning = FALSE, eval = FALSE-----------------------------
#  query_met <- GDCquery(
#      project = "TCGA-COAD",
#      data.category = "DNA Methylation",
#      platform = c("Illumina Human Methylation 450")
#  )
#  query_exp <- GDCquery(
#      project = "TCGA-COAD",
#      data.category = "Transcriptome Profiling",
#      data.type = "Gene Expression Quantification",
#      workflow.type = "STAR - Counts"
#  )
#  
#  # Get all patients that have DNA methylation and gene expression.
#  common.patients <- intersect(
#      substr(getResults(query_met, cols = "cases"), 1, 12),
#      substr(getResults(query_exp, cols = "cases"), 1, 12)
#  )
#  
#  # Only seelct the first 5 patients
#  query_met <- GDCquery(
#      project = "TCGA-COAD",
#      data.category = "DNA Methylation",
#      platform = c("Illumina Human Methylation 450"),
#      barcode = common.patients[1:5]
#  )
#  
#  query_exp <- GDCquery(
#      project = "TCGA-COAD",
#      data.category = "Transcriptome Profiling",
#      data.type = "Gene Expression Quantification",
#      workflow.type = "STAR - Counts",
#      barcode = common.patients[1:5]
#  )

## ----results_matched, message=FALSE, warning=FALSE, eval = FALSE--------------
#  datatable(
#      getResults(query_met, cols = c("data_type","cases")),
#      filter = 'top',
#      options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
#      rownames = FALSE
#  )
#  datatable(
#      getResults(query_exp, cols = c("data_type","cases")),
#      filter = 'top',
#      options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
#      rownames = FALSE
#  )

## ----message=FALSE, warning=FALSE---------------------------------------------
query <- GDCquery(
    project = "TCGA-ACC", 
    data.category = "Sequencing Reads",
    data.type = "Aligned Reads", 
    data.format = "bam",
    workflow.type = "STAR 2-Pass Transcriptome"
)
# Only first 10 to make render faster
datatable(
    getResults(query, rows = 1:10,cols = c("file_name","cases")), 
    filter = 'top',
    options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
    rownames = FALSE
)

query <- GDCquery(
    project = "TCGA-ACC", 
    data.category = "Sequencing Reads",
    data.type = "Aligned Reads", 
    data.format = "bam",
    workflow.type = "STAR 2-Pass Genome"
)
# Only first 10 to make render faster
datatable(
    getResults(query, rows = 1:10,cols = c("file_name","cases")), 
    filter = 'top',
    options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
    rownames = FALSE
)

query <- GDCquery(
    project = "TCGA-ACC", 
    data.category = "Sequencing Reads",
    data.type = "Aligned Reads", 
    data.format = "bam",
    workflow.type = "STAR 2-Pass Chimeric"
)
# Only first 10 to make render faster
datatable(
    getResults(query, rows = 1:10,cols = c("file_name","cases")), 
    filter = 'top',
    options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
    rownames = FALSE
)

query <- GDCquery(
    project = "TCGA-ACC", 
    data.category = "Sequencing Reads",
    data.type = "Aligned Reads", 
    data.format = "bam",
    workflow.type = "BWA-aln"
)
# Only first 10 to make render faster
datatable(
    getResults(query, rows = 1:10,cols = c("file_name","cases")), 
    filter = 'top',
    options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
    rownames = FALSE
)

query <- GDCquery(
    project = "TCGA-ACC", 
    data.category = "Sequencing Reads",
    data.type = "Aligned Reads", 
    data.format = "bam",
    workflow.type = "BWA with Mark Duplicates and BQSR"
)
# Only first 10 to make render faster
datatable(
    getResults(query, rows = 1:10,cols = c("file_name","cases")), 
    filter = 'top',
    options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
    rownames = FALSE
)

## ----message=FALSE, warning=FALSE---------------------------------------------
getManifest(query,save = FALSE) 

## ----message=FALSE, warning=FALSE---------------------------------------------

datatable(
    getResults(TCGAbiolinks:::GDCquery_ATAC_seq())[,c("file_name","file_size")], 
    filter = 'top',
    options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
    rownames = FALSE
)

## ----message=FALSE, warning=FALSE,eval = FALSE--------------------------------
#  query <- TCGAbiolinks:::GDCquery_ATAC_seq(file.type = "rds")
#  GDCdownload(query, method = "client")
#  
#  query <- TCGAbiolinks:::GDCquery_ATAC_seq(file.type = "bigWigs")
#  GDCdownload(query, method = "client")
#  

## ----message=FALSE, warning=FALSE,eval = TRUE---------------------------------
tab <-  getSampleFilesSummary(project = "TCGA-ACC")
datatable(
    head(tab),
    filter = 'top',
    options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
    rownames = FALSE
)

