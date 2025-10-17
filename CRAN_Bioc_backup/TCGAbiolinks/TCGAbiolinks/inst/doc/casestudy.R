## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(dpi = 300)
knitr::opts_chunk$set(cache=FALSE)

## ----message = FALSE, warning = FALSE, include = FALSE------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE-------------------------
#  library(SummarizedExperiment)
#  library(TCGAbiolinks)
#  
#  query.exp <- GDCquery(
#      project = "TCGA-BRCA",
#      data.category = "Transcriptome Profiling",
#      data.type = "Gene Expression Quantification",
#      workflow.type = "STAR - Counts",
#      sample.type = c("Primary Tumor","Solid Tissue Normal")
#  )
#  GDCdownload(
#      query = query.exp,
#      files.per.chunk = 100
#  )
#  
#  brca.exp <- GDCprepare(
#      query = query.exp,
#      save = TRUE,
#      save.filename = "brcaExp.rda"
#  )
#  
#  # get subtype information
#  infomation.subtype <- TCGAquery_subtype(tumor = "BRCA")
#  
#  # get clinical data
#  information.clinical <- GDCquery_clinic(project = "TCGA-BRCA",type = "clinical")
#  
#  # Which samples are Primary Tumor
#  samples.primary.tumour <- brca.exp$barcode[brca.exp$shortLetterCode == "TP"]
#  
#  # which samples are solid tissue normal
#  samples.solid.tissue.normal <- brca.exp$barcode[brca.exp$shortLetterCode == "NT"]

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE-------------------------
#  dataPrep <- TCGAanalyze_Preprocessing(
#      object = brca.exp,
#      cor.cut = 0.6
#  )
#  
#  dataNorm <- TCGAanalyze_Normalization(
#      tabDF = dataPrep,
#      geneInfo = geneInfoHT,
#      method = "gcContent"
#  )
#  
#  dataFilt <- TCGAanalyze_Filtering(
#      tabDF = dataNorm,
#      method = "quantile",
#      qnt.cut =  0.25
#  )
#  
#  dataDEGs <- TCGAanalyze_DEA(
#      mat1 = dataFilt[,samples.solid.tissue.normal],
#      mat2 = dataFilt[,samples.primary.tumour],
#      Cond1type = "Normal",
#      Cond2type = "Tumor",
#      fdr.cut = 0.01 ,
#      logFC.cut = 2,
#      method = "glmLRT",
#      pipeline = "edgeR"
#  )

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE-------------------------
#  ansEA <- TCGAanalyze_EAcomplete(
#      TFname = "DEA genes Normal Vs Tumor",
#      RegulonList = dataDEGs$gene_name
#  )
#  
#  TCGAvisualize_EAbarplot(
#      tf = rownames(ansEA$ResBP),
#      GOBPTab = ansEA$ResBP,
#      GOCCTab = ansEA$ResCC,
#      GOMFTab = ansEA$ResMF,
#      PathTab = ansEA$ResPat,
#      nRGTab = dataDEGs$gene_name,
#      nBar = 10
#  )

## ----fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(png)
library(grid)
img <- readPNG("case1_EA.png")
grid.raster(img)

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE-------------------------
#  
#  group1 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("NT"))
#  group2 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("TP"))
#  
#  dataSurv <- TCGAanalyze_SurvivalKM(
#      clinical_patient = dataClin,
#      dataGE = dataFilt,
#      Genelist = rownames(dataDEGs),
#      Survresult = FALSE,
#      ThreshTop = 0.67,
#      ThreshDown = 0.33,
#      p.cut = 0.05,
#      group1 = group1,
#      group2 = group2
#  )

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE-------------------------
#  
#  require(dnet)  # to change
#  org.Hs.string <- dRDataLoader(RData = "org.Hs.string")
#  
#  TabCoxNet <- TCGAvisualize_SurvivalCoxNET(
#      dataClin,
#      dataFilt,
#      Genelist = rownames(dataSurv),
#      scoreConfidence = 700,
#      org.Hs.string = org.Hs.string,
#      titlePlot = "Case Study n.1 dnet"
#  )

## ----fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(png)
library(grid)
img <- readPNG("case1_dnet.png")
grid.raster(img)

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE-------------------------
#  library(TCGAbiolinks)
#  library(SummarizedExperiment)
#  
#  query.exp <- GDCquery(
#      project = "TCGA-LGG",
#      data.category = "Transcriptome Profiling",
#      data.type = "Gene Expression Quantification",
#      workflow.type = "STAR - Counts",
#      sample.type = c("Primary Tumor")
#  )
#  
#  GDCdownload(query.exp)
#  
#  lgg.exp <- GDCprepare(
#      query = query.exp,
#      save = FALSE
#  )

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE-------------------------
#  library(dplyr)
#  
#  dataPrep <- TCGAanalyze_Preprocessing(
#      object = lgg.exp,
#      cor.cut = 0.6
#  )
#  dataNorm <- TCGAanalyze_Normalization(
#      tabDF = dataPrep,
#      geneInfo = geneInfoHT,
#      method = "gcContent"
#  )
#  
#  datFilt <- dataNorm %>%
#      TCGAanalyze_Filtering(method = "varFilter") %>%
#      TCGAanalyze_Filtering(method = "filter1") %>%
#      TCGAanalyze_Filtering(method = "filter2",foldChange = 1)
#  
#  data_Hc2 <- TCGAanalyze_Clustering(
#      tabDF = datFilt,
#      method = "consensus",
#      methodHC = "ward.D2"
#  )
#  # Add  cluster information to Summarized Experiment
#  colData(lgg.exp)$groupsHC <- paste0("EC",data_Hc2[[4]]$consensusClass)

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE-------------------------
#  TCGAanalyze_survival(
#      data = colData(lgg.exp),
#      clusterCol = "groupsHC",
#      main = "TCGA kaplan meier survival plot from consensus cluster",
#      legend = "RNA Group",
#      height = 10,
#      risk.table = T,
#      conf.int = F,
#      color = c("black","red","blue","green3"),
#      filename = "survival_lgg_expression_subtypes.png"
#  )

## ----fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(png)
library(grid)
img <- readPNG("case2_surv.png")
grid.raster(img)

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE-------------------------
#  TCGAvisualize_Heatmap(
#      data = t(datFilt),
#      col.metadata =  colData(lgg.exp)[,
#                                       c("barcode",
#                                         "groupsHC",
#                                         "paper_Histology",
#                                         "paper_IDH.codel.subtype")
#      ],
#      col.colors =  list(
#          groupsHC = c(
#              "EC1"="black",
#              "EC2"="red",
#              "EC3"="blue",
#              "EC4"="green3")
#      ),
#      sortCol = "groupsHC",
#      type = "expression", # sets default color
#      scale = "row", # use z-scores for better visualization. Center gene expression level around 0.
#      title = "Heatmap from concensus cluster",
#      filename = "case2_Heatmap.png",
#      extremes = seq(-2,2,1),
#      color.levels = colorRampPalette(c("green", "black", "red"))(n = 5),
#      cluster_rows = TRUE,
#      cluster_columns = FALSE,
#      width = 1000,
#      height = 500
#  )

## ----fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(jpeg)
library(grid)
img <- readJPEG("case2_Heatmap.jpg")
grid.raster(img)

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE-------------------------
#  library(maftools)
#  library(dplyr)
#  query <- GDCquery(
#      project = "TCGA-LGG",
#      data.category = "Simple Nucleotide Variation",
#      access = "open",
#      data.type = "Masked Somatic Mutation",
#      workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
#  )
#  GDCdownload(query)
#  LGGmut <- GDCprepare(query)
#  # Selecting gene
#  LGGmut.atr <- LGGmut %>% dplyr::filter(Hugo_Symbol == "ATR")
#  
#  dataMut <- LGGmut.atr[!duplicated(LGGmut.atr$Tumor_Sample_Barcode),]
#  dataMut$Tumor_Sample_Barcode <- substr(dataMut$Tumor_Sample_Barcode,1,12)
#  
#  # Adding the Expression Cluster classification found before
#  dataMut <- merge(dataMut, cluster, by.y = "patient", by.x = "Tumor_Sample_Barcode")

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE-------------------------
#  library(TCGAbiolinks)
#  library(SummarizedExperiment)
#  
#  #-----------------------------------
#  # STEP 1: Search, download, prepare |
#  #-----------------------------------
#  # 1.1 - DNA methylation
#  # ----------------------------------
#  query.met <- GDCquery(
#      project = "TCGA-ACC",
#      data.category = "DNA Methylation",
#      data.type = "Methylation Beta Value",
#      platform = "Illumina Human Methylation 450"
#  )
#  
#  GDCdownload(
#      query = query.met,
#      files.per.chunk = 20,
#      directory = "case3/GDCdata"
#  )
#  
#  acc.met <- GDCprepare(
#      query = query.met,
#      save = FALSE,
#      directory = "case3/GDCdata"
#  )
#  
#  #-----------------------------------
#  # 1.2 - RNA expression
#  # ----------------------------------
#  query.exp <- GDCquery(
#      project = "TCGA-ACC",
#      data.category = "Transcriptome Profiling",
#      data.type = "Gene Expression Quantification",
#      workflow.type = "STAR - Counts"
#  )
#  
#  GDCdownload(
#      query = query.exp,
#      files.per.chunk = 20,
#      directory = "case3/GDCdata"
#  )
#  
#  acc.exp <- GDCprepare(
#      query = query.exp,
#      save = FALSE,
#      directory = "case3/GDCdata"
#  )

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE-------------------------
#  # na.omit
#  acc.met <- acc.met[rowSums(is.na(assay(acc.met))) == 0,]
#  
#  # Volcano plot
#  acc.met <- TCGAanalyze_DMC(
#      data = acc.met,
#      groupCol = "subtype_MethyLevel",
#      group1 = "CIMP-high",
#      group2="CIMP-low",
#      p.cut = 10^-5,
#      diffmean.cut = 0.25,
#      legend = "State",
#      plot.filename = "case3/CIMP-highvsCIMP-low_metvolcano.png"
#  )

## ----fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(png)
library(grid)
img <- readPNG("CIMP-highvsCIMP-low_metvolcano.png")
grid.raster(img)

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE-------------------------
#  #-------------------------------------------------
#  # 2.3 - DEA - Expression analysis - volcano plot
#  # ------------------------------------------------
#  acc.exp.aux <- subset(
#      acc.exp,
#      select = colData(acc.exp)$subtype_MethyLevel %in% c("CIMP-high","CIMP-low")
#  )
#  
#  idx <- colData(acc.exp.aux)$subtype_MethyLevel %in% c("CIMP-high")
#  idx2 <- colData(acc.exp.aux)$subtype_MethyLevel %in% c("CIMP-low")
#  
#  dataPrep <- TCGAanalyze_Preprocessing(
#      object = acc.exp.aux,
#      cor.cut = 0.6
#  )
#  
#  dataNorm <- TCGAanalyze_Normalization(
#      tabDF = dataPrep,
#      geneInfo = geneInfoHT,
#      method = "gcContent"
#  )
#  
#  dataFilt <- TCGAanalyze_Filtering(
#      tabDF = dataNorm,
#      qnt.cut = 0.25,
#      method = 'quantile'
#  )
#  
#  dataDEGs <- TCGAanalyze_DEA(
#      mat1 = dataFilt[,idx],
#      mat2 = dataFilt[,idx2],
#      Cond1type = "CIMP-high",
#      Cond2type = "CIMP-low",
#      method = "glmLRT"
#  )
#  
#  TCGAVisualize_volcano(
#      x = dataDEGs$logFC,
#      y = dataDEGs$FDR,
#      filename = "case3/Case3_volcanoexp.png",
#      x.cut = 3,
#      y.cut = 10^-5,
#      names = rownames(dataDEGs),
#      color = c("black","red","darkgreen"),
#      names.size = 2,
#      xlab = " Gene expression fold change (Log2)",
#      legend = "State",
#      title = "Volcano plot (CIMP-high vs CIMP-low)",
#      width = 10
#  )

## ----fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(png)
library(grid)
img <- readPNG("figure5exp.png")
grid.raster(img)

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE-------------------------
#  #------------------------------------------
#  # 2.4 - Starburst plot
#  # -----------------------------------------
#  # If true the argument names of the genes in circle
#  # (biologically significant genes, has a change in gene
#  # expression and DNA methylation and respects all the thresholds)
#  # will be shown
#  # these genes are returned by the function see starburst object after the function is executed
#  starburst <- TCGAvisualize_starburst(
#      met = acc.met,
#      exp = dataDEGs,
#      genome = "hg19"
#      group1 = "CIMP-high",
#      group2 = "CIMP-low",
#      filename = "case3/starburst.png",
#      met.platform = "450K",
#      met.p.cut = 10^-5,
#      exp.p.cut = 10^-5,
#      diffmean.cut = 0.25,
#      logFC.cut = 3,
#      names = FALSE,
#      height = 10,
#      width = 15,
#      dpi = 300
#  )

## ----fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(png)
library(grid)
img <- readPNG("figure5star.png")
grid.raster(img)

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE-------------------------
#  library(TCGAbiolinks)
#  library(SummarizedExperiment)
#  library(ELMER)
#  library(parallel)
#  dir.create("case4")
#  setwd("case4")
#  #-----------------------------------
#  # STEP 1: Search, download, prepare |
#  #-----------------------------------
#  # 1.1 - DNA methylation
#  # ----------------------------------
#  query.met <- GDCquery(
#      project = "TCGA-KIRC",
#      data.category = "DNA Methylation",
#      data.type = "Methylation Beta Value",
#      platform = "Illumina Human Methylation 450"
#  )
#  GDCdownload(query.met)
#  kirc.met <- GDCprepare(
#      query = query.met,
#      save = TRUE,
#      save.filename = "kircDNAmet.rda",
#      summarizedExperiment = TRUE
#  )

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE-------------------------
#  # Step 1.2 download expression data
#  #-----------------------------------
#  # 1.2 - RNA expression
#  # ----------------------------------
#  query.exp <- GDCquery(
#      project = "TCGA-KIRC",
#      data.category = "Transcriptome Profiling",
#      data.type = "Gene Expression Quantification",
#      workflow.type = "STAR - Counts"
#  )
#  GDCdownload(query.exp,files.per.chunk = 20)
#  kirc.exp <- GDCprepare(
#      query = query.exp,
#      save = TRUE,
#      save.filename = "kircExp.rda"
#  )

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE-------------------------
#  distal.probes <- get.feature.probe(genome = "hg38", met.platform = "450K")

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE-------------------------
#  library(MultiAssayExperiment)
#  mae <- createMAE(
#      exp = kirc.exp,
#      met = kirc.met,
#      save = FALSE,
#      linearize.exp = TRUE,
#      filter.probes = distal.probes,
#      save.filename = "mae_kirc.rda",
#      met.platform = "450K",
#      genome = "hg38",
#      TCGA = TRUE
#  )
#  # Remove FFPE samples
#  mae <- mae[,!mae$is_ffpe]

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE-------------------------
#  group.col <- "definition"
#  group1 <-  "Primary Tumor"
#  group2 <- "Solid Tissue Normal"
#  direction <- "hypo"
#  dir.out <- file.path("kirc",direction)
#  dir.create(dir.out, recursive = TRUE)
#  #--------------------------------------
#  # STEP 3: Analysis                     |
#  #--------------------------------------
#  # Step 3.1: Get diff methylated probes |
#  #--------------------------------------
#  sig.diff <- get.diff.meth(
#      data = mae,
#      group.col = group.col,
#      group1 =  group1,
#      group2 = group2,
#      minSubgroupFrac = 0.2,
#      sig.dif = 0.3,
#      diff.dir = direction, # Search for hypomethylated probes in group 1
#      cores = 1,
#      dir.out = dir.out,
#      pvalue = 0.01
#  )
#  
#  #-------------------------------------------------------------
#  # Step 3.2: Identify significant probe-gene pairs            |
#  #-------------------------------------------------------------
#  # Collect nearby 20 genes for Sig.probes
#  nearGenes <- GetNearGenes(
#      data = mae,
#      probes = sig.diff$probe,
#      numFlankingGenes = 20, # 10 upstream and 10 dowstream genes
#      cores = 1
#  )
#  
#  pair <- get.pair(
#      data = mae,
#      group.col = group.col,
#      group1 =  group1,
#      group2 = group2,
#      nearGenes = nearGenes,
#      minSubgroupFrac = 0.4, # % of samples to use in to create groups U/M
#      permu.dir = file.path(dir.out,"permu"),
#      permu.size = 100, # Please set to 100000 to get significant results
#      raw.pvalue  = 0.05,
#      Pe = 0.01, # Please set to 0.001 to get significant results
#      filter.probes = TRUE, # See preAssociationProbeFiltering function
#      filter.percentage = 0.05,
#      filter.portion = 0.3,
#      dir.out = dir.out,
#      cores = 1,
#      label = direction
#  )
#  
#  # Identify enriched motif for significantly hypomethylated probes which
#  # have putative target genes.
#  enriched.motif <- get.enriched.motif(
#      data = mae,
#      probes = pair$Probe,
#      dir.out = dir.out,
#      label = direction,
#      min.incidence = 10,
#      lower.OR = 1.1
#  )
#  
#  TF <- get.TFs(
#      data = mae,
#      group.col = group.col,
#      group1 =  group1,
#      group2 = group2,
#      minSubgroupFrac = 0.4,
#      enriched.motif = enriched.motif,
#      dir.out = dir.out,
#      cores = 1,
#      label = direction
#  )

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE-------------------------
#  scatter.plot(
#      data = mae,
#      byProbe = list(probe = sig.diff$probe[1], numFlankingGenes = 20),
#      category = "definition",
#      dir.out = "plots",
#      lm = TRUE, # Draw linear regression curve
#      save = TRUE
#  )

## ----fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(png)
library(grid)
img <- readPNG("case4_elmer.png")
grid.raster(img)

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE-------------------------
#  scatter.plot(
#      data = mae,
#      byTF = list(
#          TF = c("RUNX1","RUNX2","RUNX3"),
#          probe = enriched.motif[[names(enriched.motif)[10]]]
#      ),
#      category = "definition",
#      dir.out = "plots",
#      save = TRUE,
#      lm_line = TRUE
#  )

## ----fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(png)
library(grid)
img <- readPNG("elmer1.png")
grid.raster(img)


## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE-------------------------
#  heatmapPairs(
#      data = mae,
#      group.col = "definition",
#      group1 = "Primary Tumor",
#      annotation.col = c("gender"),
#      group2 = "Solid Tissue Normal",
#      pairs = pair,
#      filename =  "heatmap.pdf"
#  )

## ----fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(jpeg)
library(grid)
img <- readJPEG("elmer2.jpg")
grid.raster(img)

## ----fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(png)
library(grid)
img <- readPNG("elmer3.png")
grid.raster(img)

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

