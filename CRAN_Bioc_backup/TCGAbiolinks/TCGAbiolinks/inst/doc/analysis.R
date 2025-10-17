## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(dpi = 300)
knitr::opts_chunk$set(cache=FALSE)

## ----echo = FALSE,hide=TRUE, message=FALSE,warning=FALSE----------------------
devtools::load_all(".")

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
library(SummarizedExperiment)
library(dplyr)
library(DT)

## ----eval = FALSE-------------------------------------------------------------
#  # You can define a list of samples to query and download providing relative TCGA barcodes.
#  listSamples <- c(
#      "TCGA-E9-A1NG-11A-52R-A14M-07","TCGA-BH-A1FC-11A-32R-A13Q-07",
#      "TCGA-A7-A13G-11A-51R-A13Q-07","TCGA-BH-A0DK-11A-13R-A089-07",
#      "TCGA-E9-A1RH-11A-34R-A169-07","TCGA-BH-A0AU-01A-11R-A12P-07",
#      "TCGA-C8-A1HJ-01A-11R-A13Q-07","TCGA-A7-A13D-01A-13R-A12P-07",
#      "TCGA-A2-A0CV-01A-31R-A115-07","TCGA-AQ-A0Y5-01A-11R-A14M-07"
#  )
#  
#  # Query platform Illumina HiSeq with a list of barcode
#  query <- GDCquery(
#      project = "TCGA-BRCA",
#      data.category = "Transcriptome Profiling",
#      data.type = "Gene Expression Quantification",
#      barcode = listSamples
#  )
#  
#  # Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
#  GDCdownload(query)
#  
#  # Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
#  # rsem.genes.results as values
#  BRCA.Rnaseq.SE <- GDCprepare(query)
#  
#  BRCAMatrix <- assay(BRCA.Rnaseq.SE,"unstranded")
#  # For gene expression if you need to see a boxplot correlation and AAIC plot to define outliers you can run
#  BRCA.RNAseq_CorOutliers <- TCGAanalyze_Preprocessing(BRCA.Rnaseq.SE)

## ----eval = TRUE, echo = FALSE,size = 8---------------------------------------
library(TCGAbiolinks)
dataGE <- dataBRCA[sample(rownames(dataBRCA),10),sample(colnames(dataBRCA),7)]

knitr::kable(
    dataGE[1:10,2:3], digits = 2, 
    caption = "Example of a matrix of gene expression (10 genes in rows and 2 samples in columns)",
    row.names = TRUE
)

## ----fig.width=6, fig.height=4, echo=FALSE, fig.align="center"----------------
library(png)
library(grid)
img <- readPNG("PreprocessingOutput.png")
grid.raster(img)

## ----eval = FALSE-------------------------------------------------------------
#  library(TCGAbiolinks)
#  
#  # normalization of genes
#  dataNorm <- TCGAanalyze_Normalization(
#      tabDF = BRCA.RNAseq_CorOutliers,
#      geneInfo =  geneInfoHT
#  )
#  
#  # quantile filter of genes
#  dataFilt <- TCGAanalyze_Filtering(
#      tabDF = dataNorm,
#      method = "quantile",
#      qnt.cut =  0.25
#  )
#  
#  # selection of normal samples "NT"
#  samplesNT <- TCGAquery_SampleTypes(
#      barcode = colnames(dataFilt),
#      typesample = c("NT")
#  )
#  
#  # selection of tumor samples "TP"
#  samplesTP <- TCGAquery_SampleTypes(
#      barcode = colnames(dataFilt),
#      typesample = c("TP")
#  )
#  
#  # Diff.expr.analysis (DEA)
#  dataDEGs <- TCGAanalyze_DEA(
#      mat1 = dataFilt[,samplesNT],
#      mat2 = dataFilt[,samplesTP],
#      Cond1type = "Normal",
#      Cond2type = "Tumor",
#      fdr.cut = 0.01 ,
#      logFC.cut = 1,
#      method = "glmLRT"
#  )
#  
#  # DEGs table with expression values in normal and tumor samples
#  dataDEGsFiltLevel <- TCGAanalyze_LevelTab(
#      FC_FDR_table_mRNA = dataDEGs,
#      typeCond1 = "Tumor",
#      typeCond2 = "Normal",
#      TableCond1 = dataFilt[,samplesTP],
#      TableCond2 = dataFilt[,samplesNT]
#  )
#  

## ----eval = TRUE, echo = FALSE------------------------------------------------
library(TCGAbiolinks)
dataDEGsFiltLevel$FDR <- format(dataDEGsFiltLevel$FDR, scientific = TRUE)
knitr::kable(
    dataDEGsFiltLevel[1:10,], digits = 2,
    caption = "Table of DEGs after DEA", row.names = FALSE
)

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE-------------------------
#  
#  query <- GDCquery(
#      project = "TCGA-BRCA",
#      data.category = "Transcriptome Profiling",
#      data.type = "Gene Expression Quantification",
#      workflow.type = "STAR - Counts"
#  )
#  
#  samplesDown <- getResults(query,cols=c("cases"))
#  
#  dataSmTP <- TCGAquery_SampleTypes(
#      barcode = samplesDown,
#      typesample = "TP"
#  )
#  
#  dataSmNT <- TCGAquery_SampleTypes(
#      barcode = samplesDown,
#      typesample = "NT"
#  )
#  dataSmTP_short <- dataSmTP[1:10]
#  dataSmNT_short <- dataSmNT[1:10]
#  
#  query.selected.samples <- GDCquery(
#      project = "TCGA-BRCA",
#      data.category = "Transcriptome Profiling",
#      data.type = "Gene Expression Quantification",
#      workflow.type = "STAR - Counts",
#      barcode = c(dataSmTP_short, dataSmNT_short)
#  )
#  
#  GDCdownload(
#      query = query.selected.samples
#  )
#  
#  dataPrep <- GDCprepare(
#      query = query.selected.samples,
#      save = TRUE
#  )
#  
#  dataPrep <- TCGAanalyze_Preprocessing(
#      object = dataPrep,
#      cor.cut = 0.6,
#      datatype = "HTSeq - Counts"
#  )
#  
#  dataNorm <- TCGAanalyze_Normalization(
#      tabDF = dataPrep,
#      geneInfo = geneInfoHT,
#      method = "gcContent"
#  )
#  
#  boxplot(dataPrep, outline = FALSE)
#  
#  boxplot(dataNorm, outline = FALSE)
#  
#  dataFilt <- TCGAanalyze_Filtering(
#      tabDF = dataNorm,
#      method = "quantile",
#      qnt.cut =  0.25
#  )
#  
#  dataDEGs <- TCGAanalyze_DEA(
#      mat1 = dataFilt[,dataSmTP_short],
#      mat2 = dataFilt[,dataSmNT_short],
#      Cond1type = "Normal",
#      Cond2type = "Tumor",
#      fdr.cut = 0.01 ,
#      logFC.cut = 1,
#      method = "glmLRT"
#  )
#  

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE-------------------------
#  require(TCGAbiolinks)
#  
#  query.miRNA <- GDCquery(
#      project = "TCGA-BRCA",
#      experimental.strategy = "miRNA-Seq",
#      data.category = "Transcriptome Profiling",
#      data.type = "miRNA Expression Quantification"
#  )
#  
#  GDCdownload(query = query.miRNA)
#  
#  dataAssy.miR <- GDCprepare(
#      query = query.miRNA
#  )
#  rownames(dataAssy.miR) <- dataAssy.miR$miRNA_ID
#  
#  # using read_count's data
#  read_countData <-  colnames(dataAssy.miR)[grep("count", colnames(dataAssy.miR))]
#  dataAssy.miR <- dataAssy.miR[,read_countData]
#  colnames(dataAssy.miR) <- gsub("read_count_","", colnames(dataAssy.miR))
#  
#  dataFilt <- TCGAanalyze_Filtering(
#      tabDF = dataAssy.miR,
#      method = "quantile",
#      qnt.cut =  0.25
#  )
#  
#  dataDEGs <- TCGAanalyze_DEA(
#      mat1 = dataFilt[,dataSmNT_short.miR],
#      mat2 = dataFilt[,dataSmTP_short.miR],
#      Cond1type = "Normal",
#      Cond2type = "Tumor",
#      fdr.cut = 0.01 ,
#      logFC.cut = 1,
#      method = "glmLRT"
#  )
#  

## ----eval = FALSE-------------------------------------------------------------
#  library(TCGAbiolinks)
#  # Enrichment Analysis EA
#  # Gene Ontology (GO) and Pathway enrichment by DEGs list
#  Genelist <- rownames(dataDEGsFiltLevel)
#  
#  ansEA <- TCGAanalyze_EAcomplete(
#      TFname = "DEA genes Normal Vs Tumor",
#      RegulonList = Genelist
#  )
#  
#  # Enrichment Analysis EA (TCGAVisualize)
#  # Gene Ontology (GO) and Pathway enrichment barPlot
#  
#  TCGAvisualize_EAbarplot(
#      tf = rownames(ansEA$ResBP),
#      GOBPTab = ansEA$ResBP,
#      GOCCTab = ansEA$ResCC,
#      GOMFTab = ansEA$ResMF,
#      PathTab = ansEA$ResPat,
#      nRGTab = Genelist,
#      nBar = 10
#  )
#  

## ----fig.width=6, fig.height=4, echo=FALSE, fig.align="center"----------------
library(png)
library(grid)
img <- readPNG("EAplot.png")
grid.raster(img)

## ----eval = FALSE-------------------------------------------------------------
#  clin.gbm <- GDCquery_clinic("TCGA-GBM", "clinical")
#  TCGAanalyze_survival(
#      data = clin.gbm,
#      clusterCol = "gender",
#      main = "TCGA Set\n GBM",
#      height = 10,
#      width=10
#  )

## ----fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(png)
library(grid)
img <- readPNG("case2_surv.png")
grid.raster(img)

## ----eval = FALSE-------------------------------------------------------------
#  library(TCGAbiolinks)
#  # Survival Analysis SA
#  
#  clinical_patient_Cancer <- GDCquery_clinic("TCGA-BRCA","clinical")
#  dataBRCAcomplete <- log2(BRCA_rnaseqv2)
#  
#  tokenStop <- 1
#  
#  tabSurvKMcomplete <- NULL
#  
#  for( i in 1: round(nrow(dataBRCAcomplete)/100)){
#      message( paste( i, "of ", round(nrow(dataBRCAcomplete)/100)))
#      tokenStart <- tokenStop
#      tokenStop <- 100 * i
#      tabSurvKM <- TCGAanalyze_SurvivalKM(
#          clinical_patient_Cancer,
#          dataBRCAcomplete,
#          Genelist = rownames(dataBRCAcomplete)[tokenStart:tokenStop],
#          Survresult = F,
#          ThreshTop = 0.67,
#          ThreshDown = 0.33
#      )
#  
#      tabSurvKMcomplete <- rbind(tabSurvKMcomplete,tabSurvKM)
#  }
#  
#  tabSurvKMcomplete <- tabSurvKMcomplete[tabSurvKMcomplete$pvalue < 0.01,]
#  tabSurvKMcomplete <- tabSurvKMcomplete[order(tabSurvKMcomplete$pvalue, decreasing=F),]
#  
#  tabSurvKMcompleteDEGs <- tabSurvKMcomplete[
#      rownames(tabSurvKMcomplete) %in% dataDEGsFiltLevel$mRNA,
#  ]

## ----fig.width=6, fig.height=4, echo=FALSE, fig.align="center"----------------
tabSurvKMcompleteDEGs$pvalue <- format(tabSurvKMcompleteDEGs$pvalue, scientific = TRUE)
knitr::kable(tabSurvKMcompleteDEGs[1:5,1:4], 
             digits = 2,
             caption = "Table KM-survival genes after SA",
             row.names = TRUE)
knitr::kable(tabSurvKMcompleteDEGs[1:5,5:7], 
             digits = 2,
             row.names = TRUE)

## ----eval = FALSE-------------------------------------------------------------
#  data <- TCGAanalyze_DMC(
#      data = data,
#      groupCol = "methylation_subtype",
#      group1 = "CIMP.H",
#      group2 = "CIMP.L",
#      p.cut = 10^-5,
#      diffmean.cut = 0.25,
#      legend = "State",
#      plot.filename = "coad_CIMPHvsCIMPL_metvolcano.png"
#  )

## ----fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(png)
library(grid)
img <- readPNG("figure5met.png")
grid.raster(img)

## ----fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(jpeg)
library(grid)
img <- readJPEG("case2_Heatmap.jpg")
grid.raster(img)

## ----eval = FALSE-------------------------------------------------------------
#  # normalization of genes
#  dataNorm <- TCGAbiolinks::TCGAanalyze_Normalization(dataBRCA, geneInfo)
#  
#  # quantile filter of genes
#  dataFilt <- TCGAanalyze_Filtering(
#      tabDF = dataNorm,
#      method = "quantile",
#      qnt.cut =  0.25
#  )
#  
#  # selection of normal samples "NT"
#  group1 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("NT"))
#  # selection of normal samples "TP"
#  group2 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("TP"))
#  
#  # Principal Component Analysis plot for ntop selected DEGs
#  pca <- TCGAvisualize_PCA(
#      dataFilt = dataFilt,
#      dataDEGsFiltLevel = dataDEGsFiltLevel,
#      ntopgenes = 200,
#      group1 = group1,
#      group2 =  group2
#  )

## ----fig.width=6, fig.height=4, echo=FALSE, fig.align="center"----------------
library(png)
library(grid)
img <- readPNG("PCAtop200DEGs.png")
grid.raster(img)

## ----include=FALSE,echo=FALSE, fig.height=5, message=FALSE, warning=FALSE,eval=FALSE----
#  query <- GDCquery(
#      project = "TCGA-GBM",
#      data.category = "DNA Methylation",
#      platform = "Illumina Human Methylation 27",
#      barcode = c(
#          "TCGA-02-0058-01A-01D-0186-05", "TCGA-12-1597-01B-01D-0915-05",
#          "TCGA-12-0829-01A-01D-0392-05", "TCGA-06-0155-01B-01D-0521-05",
#          "TCGA-02-0099-01A-01D-0199-05", "TCGA-19-4068-01A-01D-1228-05",
#          "TCGA-19-1788-01A-01D-0595-05", "TCGA-16-0848-01A-01D-0392-05"
#      )
#  )
#  GDCdownload(query, method = "api")
#  data <- GDCprepare(query)
#  

## ----eval=FALSE, echo=TRUE, fig.height=5, message=FALSE, warning=FALSE--------
#  query <- GDCquery(
#      project = "TCGA-GBM",
#      data.category = "DNA methylation",
#      platform = "Illumina Human Methylation 27",
#      barcode = c(
#          "TCGA-02-0058-01A-01D-0186-05", "TCGA-12-1597-01B-01D-0915-05",
#          "TCGA-12-0829-01A-01D-0392-05", "TCGA-06-0155-01B-01D-0521-05",
#          "TCGA-02-0099-01A-01D-0199-05", "TCGA-19-4068-01A-01D-1228-05",
#          "TCGA-19-1788-01A-01D-0595-05", "TCGA-16-0848-01A-01D-0392-05"
#      )
#  )
#  GDCdownload(query, method = "api")
#  data <- GDCprepare(query)

## ----eval = FALSE-------------------------------------------------------------
#  starburst <- TCGAvisualize_starburst(
#      met = coad.SummarizeExperiment,
#      exp = different.experssion.analysis.data,
#      group1 = "CIMP.H",
#      group2 = "CIMP.L",
#      met.platform = "450K",
#      genome = "hg19",
#      met.p.cut = 10^-5,
#      exp.p.cut = 10^-5,
#      names = TRUE
#  )

## ----fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(png)
library(grid)
img <- readPNG("figure5star.png")
grid.raster(img)

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

