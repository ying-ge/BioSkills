## ----label= 'Load a sample dataset', eval=TRUE--------------------------------
library(RTN)
data(tniData)

## ----label='Create a new TNI object', eval=TRUE, results='hide'---------------
# Input 1: 'expData', a named gene expression matrix (genes on rows, samples on cols); 
# Input 2: 'regulatoryElements', a vector listing genes regarded as TFs
# Input 3: 'rowAnnotation', an optional data frame with gene annotation
# Input 4: 'colAnnotation', an optional data frame with sample annotation
tfs <- c("FOXM1","E2F2","E2F3","RUNX2","PTTG1")
rtni <- tni.constructor(expData = tniData$expData, 
                        regulatoryElements = tfs, 
                        rowAnnotation = tniData$rowAnnotation, 
                        colAnnotation = tniData$colAnnotation)
# p.s. alternatively, 'expData' can be a 'SummarizedExperiment' object

## ----label='Permutation', eval=TRUE, results='hide'---------------------------
# Please set nPermutations >= 1000
rtni <- tni.permutation(rtni, nPermutations = 100)

## ----label='Bootstrap', eval=TRUE, results='hide'-----------------------------
rtni <- tni.bootstrap(rtni)

## ----label='Run DPI filter', eval=TRUE, results='hide'------------------------
rtni <- tni.dpi.filter(rtni)

## ----label='Check overall summary', eval=TRUE, collapse=TRUE------------------
tni.regulon.summary(rtni)

## ----label='Check summary of a regulon', eval=TRUE, collapse=TRUE-------------
tni.regulon.summary(rtni, regulatoryElements = "FOXM1")

## ----label='Check results', eval=TRUE, collapse=TRUE--------------------------
regulons <- tni.get(rtni, what = "regulons.and.mode", idkey = "SYMBOL")
head(regulons$FOXM1)

## ----label='Get graph', eval=TRUE---------------------------------------------
g <- tni.graph(rtni, regulatoryElements = c("FOXM1","E2F2","PTTG1"))

## ----label='Plot graph', eval=FALSE-------------------------------------------
#  library(RedeR)
#  rdp <- RedPort()
#  calld(rdp)
#  addGraph(rdp, g, layout=NULL)
#  addLegend.color(rdp, g, type="edge")
#  addLegend.shape(rdp, g)
#  relax(rdp, ps = TRUE)

## ----label='Create a new TNA object (preprocess TNI-to-TNA)', eval=TRUE, results='hide'----
# Input 1: 'object', a TNI object with regulons
# Input 2: 'phenotype', a named numeric vector, usually log2 differential expression levels
# Input 3: 'hits', a character vector, usually a set of differentially expressed genes
# Input 4: 'phenoIDs', an optional data frame with gene anottation mapped to the phenotype
data(tnaData)
rtna <- tni2tna.preprocess(object = rtni, 
                           phenotype = tnaData$phenotype, 
                           hits = tnaData$hits, 
                           phenoIDs = tnaData$phenoIDs)

## ----label='Run the MRA method', eval=TRUE, results='hide'--------------------
# Run the MRA method
rtna <- tna.mra(rtna)

## ----label='Get MRA results', eval=TRUE, collapse=TRUE------------------------
# Get MRA results;
#..setting 'ntop = -1' will return all results, regardless of a threshold
mra <- tna.get(rtna, what="mra", ntop = -1)
head(mra)

## ----label='Run GSEA method', eval=TRUE, results='hide'-----------------------
# Run the GSEA method
# Please set nPermutations >= 1000
rtna <- tna.gsea1(rtna, nPermutations=100)

## ----label='Get GSEA results', eval=TRUE, collapse=TRUE-----------------------
# Get GSEA results
gsea1 <- tna.get(rtna, what="gsea1", ntop = -1)
head(gsea1)

## ----label='Plot GSEA results', eval=FALSE------------------------------------
#  # Plot GSEA results
#  tna.plot.gsea1(rtna, labPheno="abs(log2 fold changes)", ntop = -1)

## ----label='Run the GSEA-2T method', eval=TRUE, results='hide'----------------
# Run the GSEA-2T method
# Please set nPermutations >= 1000
rtna <- tna.gsea2(rtna, nPermutations = 100)

## ----label='Get GSEA-2T results', eval=TRUE, collapse=TRUE--------------------
# Get GSEA-2T results
gsea2 <- tna.get(rtna, what = "gsea2", ntop = -1)
head(gsea2$differential)

## ----label='Plot GSEA-2T results', eval=FALSE---------------------------------
#  # Plot GSEA-2T results
#  tna.plot.gsea2(rtna, labPheno="log2 fold changes", tfs="PTTG1")

## ----eval=FALSE---------------------------------------------------------------
#  library(RTN)
#  library(TCGAbiolinks)
#  library(SummarizedExperiment)
#  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#  library(snow)

## ----eval=FALSE---------------------------------------------------------------
#  # Set GDCquery for the TCGA-BRCA cohort
#  # Gene expression data will be aligned against hg38
#  query <- GDCquery(project = "TCGA-BRCA",
#                    data.category = "Transcriptome Profiling",
#                    data.type = "Gene Expression Quantification",
#                    workflow.type = "HTSeq - FPKM-UQ",
#                    sample.type = c("Primary solid Tumor"))
#  
#  # Get a subset for demonstration (n = 500 cases)
#  cases <- getResults(query, cols = "cases")
#  cases <- sample(cases, size = 500)
#  query <- GDCquery(project = "TCGA-BRCA",
#                    data.category = "Transcriptome Profiling",
#                    data.type = "Gene Expression Quantification",
#                    workflow.type = "HTSeq - FPKM-UQ",
#                    sample.type = c("Primary solid Tumor"),
#                    barcode = cases)
#  GDCdownload(query)
#  tcgaBRCA_mRNA_data <- GDCprepare(query)

## ----eval=FALSE---------------------------------------------------------------
#  # Subset by known gene locations
#  geneRanges <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
#  tcgaBRCA_mRNA_data <- subsetByOverlaps(tcgaBRCA_mRNA_data, geneRanges)

## ----eval=FALSE---------------------------------------------------------------
#  nrow(rowData(tcgaBRCA_mRNA_data))
#  ## [1] ~30,000

## ----eval=FALSE---------------------------------------------------------------
#  # Change column names in gene annotation for better summarizations
#  colnames(rowData(tcgaBRCA_mRNA_data)) <- c("ENSEMBL", "SYMBOL", "OG_ENSEMBL")

## ----eval=FALSE---------------------------------------------------------------
#  # Save the preprocessed data for subsequent analyses
#  save(tcgaBRCA_mRNA_data, file = "tcgaBRCA_mRNA_data_preprocessed.RData")

## ----eval=FALSE---------------------------------------------------------------
#  # Load TF annotation
#  data("tfsData")

## ----eval=FALSE---------------------------------------------------------------
#  # Check TF annotation:
#  # Intersect TFs from Lambert et al. (2018) with gene annotation
#  # from the TCGA-BRCA cohort
#  geneannot <- rowData(tcgaBRCA_mRNA_data)
#  regulatoryElements <- intersect(tfsData$Lambert2018$SYMBOL, geneannot$SYMBOL)

## ----eval=FALSE---------------------------------------------------------------
#  # Run the TNI constructor
#  rtni_tcgaBRCA <- tni.constructor(expData = tcgaBRCA_mRNA_data,
#                                   regulatoryElements = regulatoryElements)

## ----eval=FALSE---------------------------------------------------------------
#  # Compute the reference regulatory network by permutation and bootstrap analyses.
#  # Please set 'spec' according to your available hardware
#  options(cluster=snow::makeCluster(spec=4, "SOCK"))
#  rtni_tcgaBRCA <- tni.permutation(rtni_tcgaBRCA, pValueCutoff = 1e-7)
#  rtni_tcgaBRCA <- tni.bootstrap(rtni_tcgaBRCA)
#  stopCluster(getOption("cluster"))

## ----eval=FALSE---------------------------------------------------------------
#  # Compute the DPI-filtered regulatory network
#  rtni_tcgaBRCA <- tni.dpi.filter(rtni_tcgaBRCA, eps = 0)

## ----eval=FALSE---------------------------------------------------------------
#  # Save the TNI object for subsequent analyses
#  save(rtni_tcgaBRCA, file="rtni_tcgaBRCA.RData")

## ----eval=FALSE---------------------------------------------------------------
#  # For example, to estimate 'alphaB' for 'nB', given 'nA', 'alphaA', and 'betaA'
#  alphaB <- tni.alpha.adjust(nB = 100, nA = 300, alphaA = 1e-5, betaA = 0.2)
#  alphaB
#  # [1] 0.029

## ----eval=FALSE---------------------------------------------------------------
#  library(RTN)
#  library(Fletcher2013b)
#  library(pheatmap)

## ----eval=FALSE---------------------------------------------------------------
#  # Load 'rtni1st' data object, which includes regulons and expression profiles
#  data("rtni1st")

## ----eval=FALSE---------------------------------------------------------------
#  # A list of transcription factors of interest (here 36 risk-associated TFs)
#  risk.tfs <- c("AFF3", "AR", "ARNT2", "BRD8", "CBFB", "CEBPB", "E2F2", "E2F3", "ENO1", "ESR1", "FOSL1", "FOXA1", "GATA3", "GATAD2A", "LZTFL1", "MTA2", "MYB", "MZF1", "NFIB", "PPARD", "RARA", "RB1", "RUNX3", "SNAPC2", "SOX10", "SPDEF", "TBX19", "TCEAL1", "TRIM29", "XBP1", "YBX1", "YPEL3", "ZNF24", "ZNF434", "ZNF552", "ZNF587")

## ----eval=FALSE---------------------------------------------------------------
#  # Compute regulon activity for individual samples
#  rtni1st <- tni.gsea2(rtni1st, regulatoryElements = risk.tfs)
#  metabric_regact <- tni.get(rtni1st, what = "regulonActivity")

## ----eval=FALSE---------------------------------------------------------------
#  # Get sample attributes from the 'rtni1st' dataset
#  metabric_annot <- tni.get(rtni1st, "colAnnotation")
#  # Get ER+/- and PAM50 attributes for pheatmap
#  attribs <- c("LumA","LumB","Basal","Her2","Normal","ER+","ER-")
#  metabric_annot <- metabric_annot[,attribs]

## ----eval=FALSE---------------------------------------------------------------
#  # Plot regulon activity profiles
#  pheatmap(t(metabric_regact$dif),
#           main="METABRIC cohort 1 (n=977 samples)",
#           annotation_col = metabric_annot,
#           show_colnames = FALSE, annotation_legend = FALSE,
#           clustering_method = "ward.D2", fontsize_row = 6,
#           clustering_distance_rows = "correlation",
#           clustering_distance_cols = "correlation")

## ----eval=FALSE, include=FALSE------------------------------------------------
#  pheatmap(t(metabric_regact$dif),
#           main="METABRIC cohort 1 (n=977 samples)",
#           treeheight_row = 40, treeheight_col = 20,
#           annotation_col = metabric_annot,
#           show_colnames = FALSE, annotation_legend = FALSE,
#           clustering_method = "ward.D2",
#           clustering_distance_rows = "correlation",
#           clustering_distance_cols = "correlation",
#           fontsize_row = 6, height = 5, width = 9,
#           filename = "fig4.pdf")

## ----eval=FALSE---------------------------------------------------------------
#  # Replace samples
#  rtni1st_tcgasamples <- tni.replace.samples(rtni1st, tcgaBRCA_mRNA_data)

## ----eval=FALSE---------------------------------------------------------------
#  # Compute regulon activity for the new samples
#  rtni1st_tcgasamples <- tni.gsea2(rtni1st_tcgasamples, regulatoryElements = risk.tfs)
#  tcga_regact <- tni.get(rtni1st_tcgasamples, what = "regulonActivity")

## ----eval=FALSE---------------------------------------------------------------
#  # Get sample attributes from the 'rtni1st_tcgasamples' dataset
#  tcga_annot <- tni.get(rtni1st_tcgasamples, "colAnnotation")
#  # Adjust PAM50 attributes for pheatmap
#  tcga_annot <- within(tcga_annot,{
#    'LumA' = ifelse(subtype_BRCA_Subtype_PAM50%in%c("LumA"),1,0)
#    'LumB' = ifelse(subtype_BRCA_Subtype_PAM50%in%c("LumB"),1,0)
#    'Basal' = ifelse(subtype_BRCA_Subtype_PAM50%in%"Basal",1,0)
#    'Her2' = ifelse(subtype_BRCA_Subtype_PAM50%in%c("Her2"),1,0)
#    'Normal' = ifelse(subtype_BRCA_Subtype_PAM50%in%c("Normal"),1,0)
#  })
#  attribs <- c("LumA","LumB","Basal","Her2","Normal")
#  tcga_annot <- tcga_annot[,attribs]

## ----eval=FALSE---------------------------------------------------------------
#  # Plot regulon activity profiles
#  pheatmap(t(tcga_regact$dif),
#           main="TCGA-BRCA cohort subset (n=500 samples)",
#           annotation_col = tcga_annot,
#           show_colnames = FALSE, annotation_legend = FALSE,
#           clustering_method = "ward.D2", fontsize_row = 6,
#           clustering_distance_rows = "correlation",
#           clustering_distance_cols = "correlation")

## ----eval=FALSE, include=FALSE------------------------------------------------
#  pheatmap(t(tcga_regact$dif),
#           main="TCGA-BRCA cohort subset (n=500 samples)",
#           treeheight_row = 40, treeheight_col = 20,
#           annotation_col = tcga_annot,
#           show_colnames = FALSE, annotation_legend = FALSE,
#           clustering_method = "ward.D2",
#           clustering_distance_rows = "correlation",
#           clustering_distance_cols = "correlation",
#           fontsize_row = 6, height = 4.5, width = 9,
#           filename = "fig5.pdf")

## ----label='Session information', eval=TRUE, echo=FALSE-----------------------
sessionInfo()

