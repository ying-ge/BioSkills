## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(dpi = 300)
knitr::opts_chunk$set(cache=FALSE)
knitr::opts_chunk$set(echo = TRUE)

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
#library(dplyr)
#library(DT)

## ----echo=TRUE, warning=FALSE, eval=TRUE--------------------------------------
projects <- getGDCprojects()$project_id
as.data.frame(grep("TARGET", projects,value = TRUE))

## ----echo=TRUE, results='hide', message=FALSE, warning=FALSE, eval=FALSE------
#  #Downloading and prepare TARGET CASE
#  query_Target <- GDCquery(project = "TARGET-AML",
#                           data.category = "Transcriptome Profiling",
#                           data.type = "Gene Expression Quantification",
#                           workflow.type = "HTSeq - Counts")
#  
#  samplesDown_Target <- getResults(query_Target, cols = c("cases"))
#  
#  dataSmTB_Target <- TCGAquery_SampleTypes(barcode = samplesDown_Target,
#                                           typesample = "TB")
#  
#  dataSmNB_Target <- TCGAquery_SampleTypes(barcode = samplesDown_Target,
#                                           typesample = "TBM")
#  dataSmTB_short_Target <- dataSmTB_Target[1:10]
#  dataSmNB_short_Target <- dataSmNB_Target[1:10]
#  
#  queryDown_Target <- GDCquery(project = "TARGET-AML",
#                               data.category = "Transcriptome Profiling",
#                               data.type = "Gene Expression Quantification",
#                               workflow.type = "HTSeq - Counts",
#                               barcode = c(dataSmTB_short_Target, dataSmNB_short_Target))
#  
#  GDCdownload(queryDown_Target)
#  
#  ### SummarizedExperiment = TRUE
#  data <- GDCprepare(queryDown_Target)
#  
#  dataPrep_Target <- TCGAanalyze_Preprocessing(object = data,
#                                               cor.cut = 0.6,
#                                               datatype = "HTSeq - Counts")
#  
#  dataNorm_Target <- TCGAanalyze_Normalization(tabDF = dataPrep_Target,
#                                               geneInfo = geneInfoHT,
#                                               method = "gcContent")
#  
#  boxplot(dataPrep_Target, outline = FALSE)
#  
#  boxplot(dataNorm_Target, outline = FALSE)
#  
#  dataFilt_Target <- TCGAanalyze_Filtering(tabDF = dataNorm_Target,
#                                           method = "quantile",
#                                           qnt.cut =  0.25)

## ----echo=TRUE, results='hide', message=FALSE, warning=FALSE, eval=FALSE------
#  CancerProject <- "TCGA-BRCA"
#  DataDirectory <- paste0("../GDC/",gsub("-","_",CancerProject))
#  FileNameData <- paste0(DataDirectory, "_","HTSeq_Counts",".rda")
#  
#  query <- GDCquery(project = CancerProject,
#                    data.category = "Transcriptome Profiling",
#                    data.type = "Gene Expression Quantification",
#                    workflow.type = "HTSeq - Counts")
#  
#  samplesDown <- getResults(query,cols = c("cases"))
#  
#  dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
#                                    typesample = "TP")
#  
#  dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
#                                    typesample = "NT")
#  dataSmTP_short <- dataSmTP[1:10]
#  dataSmNT_short <- dataSmNT[1:10]
#  
#  queryDown <- GDCquery(project = CancerProject,
#                        data.category = "Transcriptome Profiling",
#                        data.type = "Gene Expression Quantification",
#                        workflow.type = "HTSeq - Counts",
#                        barcode = c(dataSmTP_short, dataSmNT_short))
#  
#  GDCdownload(query = queryDown)
#  
#  dataPrep1 <- GDCprepare(query = queryDown,
#                          save = TRUE,
#                          save.filename = "TCGA_BRCA_HTSeq_Countds.rda")
#  
#  dataPrep <- TCGAanalyze_Preprocessing(object = dataPrep1,
#                                        cor.cut = 0.6,
#                                        datatype = "HTSeq - Counts")
#  
#  
#  dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
#                                        geneInfo = geneInfoHT,
#                                        method = "gcContent")
#  dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
#                                        geneInfo = geneInfoHT,
#                                        method = "geneLength")
#  
#  dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
#                                    method = "quantile",
#                                    qnt.cut =  0.25)

## ----echo=TRUE, results='hide', message=FALSE, warning=FALSE, eval=FALSE------
#  ### dataframe will have genes from dataFilt but raw counts from dataPrep
#  dataPrep_raw <- UseRaw_afterFilter(dataPrep, dataFilt)

## ----echo=TRUE, results='hide', message=FALSE, warning=FALSE, eval=FALSE------
#  ##use previously fetched BRCA data:
#  Pam50.subtypes<-TCGA_MolecularSubtype(colnames(dataFilt))

## ----echo=TRUE, results='hide', message=FALSE, warning=FALSE, eval=FALSE------
#  dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,dataSmTP_short],
#                              mat2 = dataFilt[,dataSmNT_short],
#                              pipeline="limma",
#                              Cond1type = "Normal",
#                              Cond2type = "Tumor",
#                              fdr.cut = 0.01 ,
#                              logFC.cut = 1,
#                              method = "glmLRT", ClinicalDF = data.frame())

## ----echo=TRUE, results='hide', message=FALSE, warning=FALSE, eval=FALSE------
#  ### ##download and prepare lung data through GDC### ##
#  query.lung <- GDCquery(project = "TCGA-LUSC",
#                         data.category = "Transcriptome Profiling",
#                         data.type = "Gene Expression Quantification",
#                         workflow.type = "HTSeq - Counts")
#  
#  samplesDown.lusc <- getResults(query.lung,cols=c("cases"))
#  
#  dataSmTP.lusc <- TCGAquery_SampleTypes(barcode = samplesDown.lusc,
#                                         typesample = "TP")
#  
#  dataSmNT.lusc <- TCGAquery_SampleTypes(barcode = samplesDown.lusc,
#                                         typesample = "NT")
#  dataSmTP_short.lusc <- dataSmTP.lusc[1:10]
#  dataSmNT_short.lusc <- dataSmNT.lusc[1:10]
#  length(dataSmTP.lusc)
#  
#  
#  queryDown.lung <- GDCquery(project = "TCGA-LUSC",
#                             data.category = "Transcriptome Profiling",
#                             data.type = "Gene Expression Quantification",
#                             workflow.type = "HTSeq - Counts",
#                             barcode = c(dataSmTP.lusc, dataSmNT.lusc))
#  
#  GDCdownload(query = queryDown.lung)
#  
#  dataPrep1.tcga <- GDCprepare(query = queryDown.lung,
#                               save = TRUE,
#                               save.filename = "TCGA_LUSC_HTSeq_Countds.rda")
#  dataPrep.tcga <- TCGAanalyze_Preprocessing(object = dataPrep1.tcga,
#                                             cor.cut = 0.6,
#                                             datatype = "HTSeq - Counts")
#  
#  dataNorm.tcga <- TCGAanalyze_Normalization(tabDF = dataPrep.tcga,
#                                             geneInfo = geneInfoHT,
#                                             method = "gcContent")
#  
#  dataNorm.tcga <- TCGAanalyze_Normalization(tabDF = dataNorm.tcga,
#                                             geneInfo = geneInfoHT,
#                                             method = "geneLength")
#  dataFilt.tcga <- TCGAanalyze_Filtering(tabDF = dataNorm.tcga,
#                                         method = "quantile",
#                                         qnt.cut =  0.25)
#  
#  ### #Filtering data so all samples have a pam50 subtype for LUSC
#  
#  diff<-setdiff(colnames(dataFilt.tcga), TCGA_MolecularSubtype(colnames(dataFilt.tcga))$filtered)
#  TCGA_MolecularSubtype(diff)$subtypes$barcodes
#  
#  dataFilt.tcga.lusc<-dataFilt.tcga[,diff]
#  LUSC.subtypes<-TCGA_MolecularSubtype(colnames(dataFilt.tcga.lusc))$subtypes$subtype
#  
#  ### Differential expression analysis after correcting for "Plate" factor.
#  DEG.lusc <- TCGAanalyze_DEA(MAT=dataFilt.tcga.lusc,
#                              pipeline="limma",
#                              batch.factors = c("Plate"),
#                              Cond1type = "Normal",
#                              Cond2type = "Tumor",
#                              fdr.cut = 0.01 ,
#                              logFC.cut = 1,
#                              method = "glmLRT", ClinicalDF = data.frame(),
#                              Condtypes = LUSC.subtypes,
#                              contrast.formula = "LUSC.basalvsLUSC.classical=LUSC.basal - LUSC.classical")

## ----echo=TRUE, results='hide', message=FALSE, warning=FALSE, eval=FALSE------
#  
#  batch.corrected <- TCGAbatch_Correction(tabDF = dataFilt.tcga, batch.factor = "Plate")
#  
#  batch.corrected.adjusted <- TCGAbatch_Correction(tabDF = dataFilt.tcga,
#                                                   batch.factor = "Plate",
#                                                   adjustment=c("TSS"))

## ----echo=TRUE, results='hide', message=FALSE, warning=FALSE, eval=FALSE------
#  
#  #dataFilt is a gene expression matrix from pancancer33 RNASeq data that can be generated following our previous described workflow
#  
#  sampleTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), typesample = "TP")
#  sampleNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), typesample = "NT")
#  dataFiltTP <- dataFilt[,sampleTP]
#  dataFiltNT <- dataFilt[,sampleNT]
#  
#  dataSubt_LUAD <-TCGAquery_subtype(tumor = "LUAD")
#  sampleStageI_LUAD <- dataSubt_LUAD[dataSubt_LUAD$Tumor.stage %in% c("Stage I","Stage IA","Stage IB"),]
#  dataFilt.LUAD.stageI.TP <- dataFiltTP[,substr(colnames(dataFiltTP),1,12) %in% sampleStageI_LUAD$patient]
#  dataFilt.LUAD.stageI.NT <- dataFiltNT[,substr(colnames(dataFiltNT),1,12) %in% sampleStageI_LUAD$patient]
#  
#  # here for example the end-user can replace LUSC with unpublished data in a dataframe containing counts.
#  # and adding the information for LUSC.stageI.TP or LUSC.stageI.NT with unpublished data
#  
#  dataSubt_LUSC <-TCGAquery_subtype(tumor = "LUSC")
#  sampleStageI_LUSC <- dataSubt_LUSC[dataSubt_LUSC$T.stage %in% c("T1","T1a","T1b"),]
#  dataFilt.LUSC.stageI.TP <- dataFiltTP[,substr(colnames(dataFiltTP),1,12) %in% sampleStageI_LUSC$patient]
#  dataFilt.LUSC.stageI.NT <- dataFiltNT[,substr(colnames(dataFiltNT),1,12) %in% sampleStageI_LUSC$patient]
#  
#  countsTable <-cbind(dataFilt.LUAD.stageI.TP,dataFilt.LUAD.stageI.NT,
#                      dataFilt.LUSC.stageI.TP, dataFilt.LUSC.stageI.NT)
#  
#  AnnotationCounts <- matrix(0,ncol(countsTable),2)
#  colnames(AnnotationCounts) <- c("Samples","Batch")
#  rownames(AnnotationCounts) <- colnames(countsTable)
#  
#  AnnotationCounts <- as.data.frame(AnnotationCounts)
#  AnnotationCounts$Samples <- colnames(countsTable)
#  
#  AnnotationCounts[colnames(dataFilt.LUAD.stageI.TP),"Batch"] <- "LUAD.stageI.TP"
#  AnnotationCounts[colnames(dataFilt.LUAD.stageI.NT),"Batch"] <- "LUAD.stageI.NT"
#  AnnotationCounts[colnames(dataFilt.LUSC.stageI.TP),"Batch"] <- "LUSC.stageI.TP"
#  AnnotationCounts[colnames(dataFilt.LUSC.stageI.NT),"Batch"] <- "LUSC.stageI.NT"
#  
#  countsCorrected <- TCGAbatch_Correction(tabDF =countsTable,
#                                          UnpublishedData = TRUE,
#                                          AnnotationDF = AnnotationCounts)

## ----echo=TRUE, results='hide', message=FALSE, warning=FALSE, eval=FALSE------
#  Purity.BRCA <- TCGAtumor_purity(colnames(dataPrep.tcga), 0, 0, 0, 0, 0.7)

## ----echo=TRUE, results='hide', message=FALSE, warning=FALSE, eval=FALSE------
#  ##Brain data through Recount2
#  
#  brain.rec <- TCGAquery_recount2(project = "gtex", tissue = "brain")

