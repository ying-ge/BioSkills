## ----eval=FALSE---------------------------------------------------------------
#  if (!require("BiocManager"))
#      install.packages("BiocManager")
#  BiocManager::install("maftools")

## ----loadlib, results='hide', message=FALSE-----------------------------------
library(maftools)

## ----readmaf------------------------------------------------------------------
#path to TCGA LAML MAF file
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') 
#clinical information containing survival information and histology. This is optional
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') 

laml = read.maf(maf = laml.maf, clinicalData = laml.clin)

## ----mafobject----------------------------------------------------------------
#Typing laml shows basic summary of MAF file.
laml

## ----mafsummary, eval=FALSE---------------------------------------------------
#  #Shows sample summry.
#  getSampleSummary(laml)
#  #Shows gene summary.
#  getGeneSummary(laml)
#  #shows clinical data associated with samples
#  getClinicalData(laml)
#  #Shows all fields in MAF
#  getFields(laml)
#  #Writes maf summary to an output file with basename laml.
#  write.mafSummary(maf = laml, basename = 'laml')

## ----summaryPlot,fig.height=4, fig.width=6------------------------------------
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

## ----oncoplot, fig.align='left',fig.height=3.5,fig.width=6, fig.align='left'----
#oncoplot for top ten mutated genes.
oncoplot(maf = laml, top = 10)

## ----titv, fig.height=3, fig.width=4.2, eval = T, fig.align='left'------------
laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = laml.titv)

## ----lollipopPlot,fig.align='left', fig.width=4.5, fig.height=2.5-------------
#lollipop plot for DNMT3A, which is one of the most frequent mutated gene in Leukemia.
lollipopPlot(
  maf = laml,
  gene = 'DNMT3A',
  AACol = 'Protein_Change',
  showMutationRate = TRUE,
  labelPos = 882
)

## ----lollipopPlot4, fig.align='left', fig.width=4.5, fig.height=2.5-----------
#example data
my_data = data.frame(pos = sample.int(912, 15, replace = TRUE), count = sample.int(30, 15, replace = TRUE))
head(my_data)
lollipopPlot(data = my_data, gene = "DNMT3A")

## ----plotProtein,fig.align='left', fig.width=5, fig.height=1.2----------------
plotProtein(gene = "TP53", refSeqID = "NM_000546")

## ----rainfallPlot, results='hide', message=FALSE------------------------------
brca <- system.file("extdata", "brca.maf.gz", package = "maftools")
brca = read.maf(maf = brca, verbose = FALSE)

## ----rainfallPlot2, fig.height=3,fig.width=6,fig.align='left'-----------------
rainfallPlot(maf = brca, detectChangePoints = TRUE, pointSize = 0.4)

## ----tcgaCompare, fig.align='left', fig.height=3.25, fig.width=6, message=FALSE, results='hide'----
laml.mutload = tcgaCompare(maf = laml, cohortName = 'Example-LAML', logscale = TRUE, capture_size = 50)

## ----plotVaf, fig.align='left', fig.height=3, fig.width=3---------------------
plotVaf(maf = laml, vafCol = 'i_TumorVAF_WU')

## ----readGistic---------------------------------------------------------------
gistic_res_folder <- system.file("extdata", package = "maftools")
laml.gistic = readGistic(gisticDir = gistic_res_folder, isTCGA = TRUE)

#GISTIC object
laml.gistic

## ----gisticChromPlot, fig.width=4, fig.height=3, fig.align='left'-------------
gisticChromPlot(gistic = laml.gistic, markBands = "all")

## ----coGisticChromPlot, fig.width=4, fig.height=5, fig.align='left'-----------
coGisticChromPlot(gistic1 = laml.gistic, gistic2 = laml.gistic, g1Name = "AML-1", g2Name = "AML-2", type = 'Amp')

## ----gisticBubblePlot, fig.width=4, fig.height=3, fig.align='left'------------
gisticBubblePlot(gistic = laml.gistic)

## ----gisticOncoPlot, fig.align='left',fig.width=5, fig.height=3, eval=T-------
gisticOncoPlot(gistic = laml.gistic, clinicalData = getClinicalData(x = laml), clinicalFeatures = 'FAB_classification', sortByAnnotation = TRUE, top = 10)

## ----plotCBSsegments, fig.height=2.5,fig.width=4,fig.align='left'-------------
tcga.ab.009.seg <- system.file("extdata", "TCGA.AB.3009.hg19.seg.txt", package = "maftools")
plotCBSsegments(cbsFile = tcga.ab.009.seg)

## ----somaticInteractions, message=FALSE, fig.height=5, fig.width=5------------
#exclusive/co-occurance event analysis on top 10 mutated genes. 
somaticInteractions(maf = laml, top = 25, pvalue = c(0.05, 0.1))

## ----oncodrive, fig.align='default', fig.width=7,fig.height=5, message=F,results='hide', eval=T----
laml.sig = oncodrive(maf = laml, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')

## -----------------------------------------------------------------------------
head(laml.sig)

## ----plotOncodrive, fig.align='left', fig.width=3.2, fig.height=3.2-----------
plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)

## ----pfamDomains, fig.align='left', fig.width=4, fig.height=3-----------------
laml.pfam = pfamDomains(maf = laml, AACol = 'Protein_Change', top = 10)
#Protein summary (Printing first 7 columns for display convenience)
laml.pfam$proteinSummary[,1:7, with = FALSE]
#Domain summary (Printing first 3 columns for display convenience)
laml.pfam$domainSummary[,1:3, with = FALSE]

## ----mafSurvival, fig.width=3, fig.height=3-----------------------------------
#Survival analysis based on grouping of DNMT3A mutation status
mafSurvival(maf = laml, genes = 'DNMT3A', time = 'days_to_last_followup', Status = 'Overall_Survival_Status', isTCGA = TRUE)

## ----survGroup----------------------------------------------------------------
#Using top 20 mutated genes to identify a set of genes (of size 2) to predict poor prognostic groups
prog_geneset = survGroup(maf = laml, top = 20, geneSetSize = 2, time = "days_to_last_followup", Status = "Overall_Survival_Status", verbose = FALSE)

print(prog_geneset)

## ----mafSurvGroup, fig.width=3, fig.height=3----------------------------------
mafSurvGroup(maf = laml, geneSet = c("DNMT3A", "FLT3"), time = "days_to_last_followup", Status = "Overall_Survival_Status")

## ----results='hide', message=FALSE--------------------------------------------
#Primary APL MAF
primary.apl = system.file("extdata", "APL_primary.maf.gz", package = "maftools")
primary.apl = read.maf(maf = primary.apl)
#Relapse APL MAF
relapse.apl = system.file("extdata", "APL_relapse.maf.gz", package = "maftools")
relapse.apl = read.maf(maf = relapse.apl)

## ----mafCompare, fig.align='left'---------------------------------------------
#Considering only genes which are mutated in at-least in 5 samples in one of the cohort to avoid bias due to genes mutated in single sample.
pt.vs.rt <- mafCompare(m1 = primary.apl, m2 = relapse.apl, m1Name = 'Primary', m2Name = 'Relapse', minMut = 5)
print(pt.vs.rt)

## ----forestPlot, fig.width=6, fig.height=4.5, fig.align='left'----------------
forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.1)

## ----coOncoplot, fig.height=2.5,fig.width=6, eval=T, fig.align='left'---------
genes = c("PML", "RARA", "RUNX1", "ARID1B", "FLT3")
coOncoplot(m1 = primary.apl, m2 = relapse.apl, m1Name = 'PrimaryAPL', m2Name = 'RelapseAPL', genes = genes, removeNonMutated = TRUE)

## ----coBarplot, fig.height=3, fig.width=4-------------------------------------
coBarplot(m1 = primary.apl, m2 = relapse.apl, m1Name = "Primary", m2Name = "Relapse")

## ----lollipopPlot2, warning=FALSE, message=FALSE,fig.align='left', results='hide', fig.height=3.5, fig.width=5----
lollipopPlot2(m1 = primary.apl, m2 = relapse.apl, gene = "PML", AACol1 = "amino_acid_change", AACol2 = "amino_acid_change", m1_name = "Primary", m2_name = "Relapse")

## ----clinicalEnrichment-------------------------------------------------------
fab.ce = clinicalEnrichment(maf = laml, clinicalFeature = 'FAB_classification')

#Results are returned as a list. Significant associations p-value < 0.05
fab.ce$groupwise_comparision[p_value < 0.05]

## ----plotEnrichmentResults, fig.width=4, fig.height=3-------------------------
plotEnrichmentResults(enrich_res = fab.ce, pVal = 0.05, geneFontSize = 0.5, annoFontSize = 0.6)

## ----drugInteractions, fig.height=3, fig.width=5------------------------------
dgi = drugInteractions(maf = laml, fontSize = 0.75)

## ----drugInteractions2--------------------------------------------------------
dnmt3a.dgi = drugInteractions(genes = "DNMT3A", drugs = TRUE)
#Printing selected columns.
dnmt3a.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]

## ----OncogenicPathways, fig.width=7, fig.height=6-----------------------------
pws = pathways(maf = laml, plotType = 'treemap')

## ----PlotOncogenicPathways, fig.width=6, fig.height=3.5-----------------------
plotPathways(maf = laml, pathlist = pws)

## ----inferHeterogeneity, echo = TRUE, fig.align='left', fig.height=3.5, fig.width=4, eval=T----
#Heterogeneity in sample TCGA.AB.2972
library("mclust")
tcga.ab.2972.het = inferHeterogeneity(maf = laml, tsb = 'TCGA-AB-2972', vafCol = 'i_TumorVAF_WU')
print(tcga.ab.2972.het$clusterMeans)
#Visualizing results
plotClusters(clusters = tcga.ab.2972.het)

## ----plotClusters, fig.align='left', fig.height=3.5, fig.width=4, eval=T------
seg = system.file('extdata', 'TCGA.AB.3009.hg19.seg.txt', package = 'maftools')
tcga.ab.3009.het = inferHeterogeneity(maf = laml, tsb = 'TCGA-AB-3009', segFile = seg, vafCol = 'i_TumorVAF_WU')
#Visualizing results. Highlighting those variants on copynumber altered variants.
plotClusters(clusters = tcga.ab.3009.het, genes = 'CN_altered', showCNvars = TRUE)

## ----trinucleotideMatrix, eval=TRUE-------------------------------------------
#Requires BSgenome object
library("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
laml.tnm = trinucleotideMatrix(maf = laml, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")

## ----plotApobecDiff, eval=TRUE, fig.height=3, fig.width=5---------------------
plotApobecDiff(tnm = laml.tnm, maf = laml, pVal = 0.2)

## ----echo=FALSE---------------------------------------------------------------
par(mar = c(2, 2, 2, 1))
plot(NA, xlim = c(1, 10), ylim = c(0, 30), frame.plot = FALSE, axes = FALSE, xlab = NA, ylab = NA)
rect(xleft = 3, ybottom = 28, xright = 7, ytop = 30, col = grDevices::adjustcolor("gray70", alpha.f = 0.6), lwd = 1.2, border = "maroon")
text(x = 5, y = 29, labels = "MAF", font = 2)
arrows(x0 = 5, y0 = 28, x1 = 5, y1 = 26, length = 0.1, lwd = 2)
text(x = 5, y = 25, labels = "trinucleotideMatrix()", font = 3)
arrows(x0 = 5, y0 = 24, x1 = 5, y1 = 21, length = 0.1, lwd = 2)
text(x = 5, y = 20, labels = "estimateSignatures()", font = 3)
arrows(x0 = 5, y0 = 19, x1 = 5, y1 = 16, length = 0.1, lwd = 2)
text(x = 5, y = 15, labels = "plotCophenetic()", font = 3)
arrows(x0 = 5, y0 = 14, x1 = 5, y1 = 11, length = 0.1, lwd = 2)
text(x = 5, y = 10, labels = "extractSignatures()", font = 3)
arrows(x0 = 5, y0 = 9, x1 = 5, y1 = 6, length = 0.1, lwd = 2)
text(x = 5, y = 5, labels = "compareSignatures()", font = 3)
arrows(x0 = 5, y0 = 4, x1 = 5, y1 = 1, length = 0.1, lwd = 2)
text(x = 5, y = 0, labels = "plotSignatures()", font = 3)

## ----estimateSignatures, fig.height=5, fig.width=5, eval=FALSE, message=FALSE----
#  library('NMF')
#  laml.sign = estimateSignatures(mat = laml.tnm, nTry = 6)

## ----estimateSignatures2, fig.height=3, fig.width=3, eval=TRUE, message=FALSE, echo=FALSE, include=FALSE----
#Run main function with maximum 6 signatures. 
library('NMF')
laml.sign = estimateSignatures(mat = laml.tnm, nTry = 6, pConstant = 0.1, plotBestFitRes = FALSE, parallel = 2)

## ----plotCophenetic, fig.width=3, fig.height=3, eval=TRUE---------------------
plotCophenetic(res = laml.sign)

## ----extractSignatures, eval=FALSE--------------------------------------------
#  laml.sig = extractSignatures(mat = laml.tnm, n = 3)

## ----extractSignatures2, eval=TRUE, echo=FALSE--------------------------------
laml.sig = extractSignatures(mat = laml.tnm, n = 3, pConstant = 0.1,  parallel = 2)

## ----compareSignatures, eval=TRUE---------------------------------------------
#Compate against original 30 signatures 
laml.og30.cosm = compareSignatures(nmfRes = laml.sig, sig_db = "legacy")
#Compate against updated version3 60 signatures 
laml.v3.cosm = compareSignatures(nmfRes = laml.sig, sig_db = "SBS")

## ----pheatmap, fig.width=7, fig.height=2.5, fig.align='center', eval=TRUE-----
library('pheatmap')
pheatmap::pheatmap(mat = laml.og30.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")

## ----plotSignatures, fig.width=6, fig.height=4, fig.align='center', eval = TRUE----
maftools::plotSignatures(nmfRes = laml.sig, title_size = 1.2, sig_db = "SBS")

## ----legoplot3d, eval=FALSE---------------------------------------------------
#  library("barplot3d")
#  #Visualize first signature
#  sig1 = laml.sig$signatures[,1]
#  barplot3d::legoplot3d(contextdata = sig1, labels = FALSE, scalexy = 0.01, sixcolors = "sanger", alpha = 0.5)

## ----annovarToMaf, eval=TRUE--------------------------------------------------
var.annovar = system.file("extdata", "variants.hg19_multianno.txt", package = "maftools")
var.annovar.maf = annovarToMaf(annovar = var.annovar, Center = 'CSI-NUS', refBuild = 'hg19', 
                               tsbCol = 'Tumor_Sample_Barcode', table = 'ensGene')


## ----icgcSimpleMutationToMAF--------------------------------------------------
#Read sample ICGC data for ESCA
esca.icgc <- system.file("extdata", "simple_somatic_mutation.open.ESCA-CN.sample.tsv.gz", package = "maftools")
esca.maf <- icgcSimpleMutationToMAF(icgc = esca.icgc, addHugoSymbol = TRUE)
#Printing first 16 columns for display convenience.
print(esca.maf[1:5,1:16, with = FALSE])

## ----prepareMutSig, eval=FALSE------------------------------------------------
#  laml.mutsig.corrected = prepareMutSig(maf = laml)
#  # Converting gene names for 1 variants from 1 genes
#  #    Hugo_Symbol MutSig_Synonym N
#  # 1:    ARHGAP35          GRLF1 1
#  # Original symbols are preserved under column OG_Hugo_Symbol.

## ----subsetMaf----------------------------------------------------------------
#Extract data for samples 'TCGA.AB.3009' and 'TCGA.AB.2933'  (Printing just 5 rows for display convenience)
subsetMaf(maf = laml, tsb = c('TCGA-AB-3009', 'TCGA-AB-2933'), mafObj = FALSE)[1:5]
##Same as above but return output as an MAF object (Default behaviour)
subsetMaf(maf = laml, tsb = c('TCGA-AB-3009', 'TCGA-AB-2933'))

## ----subsetMaf2---------------------------------------------------------------
#Select all Splice_Site mutations from DNMT3A and NPM1
subsetMaf(maf = laml, genes = c('DNMT3A', 'NPM1'), mafObj = FALSE,query = "Variant_Classification == 'Splice_Site'")

#Same as above but include only 'i_transcript_name' column in the output
subsetMaf(maf = laml, genes = c('DNMT3A', 'NPM1'), mafObj = FALSE, query = "Variant_Classification == 'Splice_Site'", fields = 'i_transcript_name')

## ----subsetMaf3---------------------------------------------------------------
#Select all samples with FAB clasification M4 in clinical data 
laml_m4 = subsetMaf(maf = laml, clinQuery = "FAB_classification %in% 'M4'")

## ----sampleSwaps, eval=FALSE--------------------------------------------------
#  #Path to BAM files
#  bams = c(
#    "DBW-40-N.bam",
#    "DBW-40-1T.bam",
#    "DBW-40-2T.bam",
#    "DBW-40-3T.bam",
#    "DBW-43-N.bam",
#    "DBW-43-1T.bam"
#  )
#  
#  res = maftools::sampleSwaps(bams = bams, build = "hg19")
#  # Fetching readcounts from BAM files..
#  # Summarizing allele frequncy table..
#  # Performing pairwise comparison..
#  # Done!

## ----pairwise_comparison, eval=FALSE------------------------------------------
#   res$pairwise_comparison

## -----------------------------------------------------------------------------
# X_bam     Y_bam concordant_snps discordant_snps fract_concordant_snps  cor_coef XY_possibly_paired
#  1: DBW-40-1T DBW-40-2T            5488             571             0.9057600 0.9656484                Yes
#  2: DBW-40-1T DBW-40-3T            5793             266             0.9560984 0.9758083                Yes
#  3: DBW-40-1T  DBW-43-N            5534             525             0.9133520 0.9667620                Yes
#  4: DBW-40-2T DBW-40-3T            5853             206             0.9660010 0.9817475                Yes
#  5: DBW-40-2T  DBW-43-N            5131             928             0.8468394 0.9297096                Yes
#  6: DBW-40-3T  DBW-43-N            5334             725             0.8803433 0.9550670                Yes
#  7:  DBW-40-N DBW-43-1T            5709             350             0.9422347 0.9725684                Yes
#  8: DBW-40-1T  DBW-40-N            2829            3230             0.4669087 0.3808831                 No
#  9: DBW-40-1T DBW-43-1T            2796            3263             0.4614623 0.3755364                 No
# 10: DBW-40-2T  DBW-40-N            2760            3299             0.4555207 0.3641647                 No
# 11: DBW-40-2T DBW-43-1T            2736            3323             0.4515597 0.3579747                 No
# 12: DBW-40-3T  DBW-40-N            2775            3284             0.4579964 0.3770581                 No
# 13: DBW-40-3T DBW-43-1T            2753            3306             0.4543654 0.3721022                 No
# 14:  DBW-40-N  DBW-43-N            2965            3094             0.4893547 0.3839140                 No
# 15: DBW-43-1T  DBW-43-N            2876            3183             0.4746658 0.3797829                 No

## ----BAM_matches, eval=FALSE--------------------------------------------------
#  res$BAM_matches

## -----------------------------------------------------------------------------
# [[1]]
# [1] "DBW-40-1T" "DBW-40-2T" "DBW-40-3T" "DBW-43-N" 
# 
# [[2]]
# [1] "DBW-40-2T" "DBW-40-3T" "DBW-43-N" 
# 
# [[3]]
# [1] "DBW-40-3T" "DBW-43-N" 
# 
# [[4]]
# [1] "DBW-40-N"  "DBW-43-1T"

## ----cor, eval=FALSE----------------------------------------------------------
#  cor_table = cor(res$AF_table)

## ----cortable, echo=FALSE-----------------------------------------------------
cor_table = structure(c(1, 0.971671361359591, 0.982608032160979, 0.381966662753787, 
0.380812832617918, 0.97246945978859, 0.971671361359591, 1, 0.988268712494756, 
0.365843547670381, 0.364768267525972, 0.933880555972565, 0.982608032160979, 
0.988268712494756, 1, 0.376165783421095, 0.375700923176265, 0.959717559987264, 
0.381966662753787, 0.365843547670381, 0.376165783421095, 1, 0.979259788301874, 
0.385257303482557, 0.380812832617918, 0.364768267525972, 0.375700923176265, 
0.979259788301874, 1, 0.384975386482912, 0.97246945978859, 0.933880555972565, 
0.959717559987264, 0.385257303482557, 0.384975386482912, 1), .Dim = c(6L, 
6L), .Dimnames = list(c("DBW-40-1T", "DBW-40-2T", "DBW-40-3T", 
"DBW-40-N", "DBW-43-1T", "DBW-43-N"), c("DBW-40-1T", "DBW-40-2T", 
"DBW-40-3T", "DBW-40-N", "DBW-43-1T", "DBW-43-N")))

## ----pheatmap2, fig.width=6, fig.height=4-------------------------------------
pheatmap::pheatmap(cor_table, breaks = seq(0, 1, 0.01))

## ----tcgaAvailable------------------------------------------------------------
tcga_avail = tcgaAvailable()
head(tcga_avail, 3)

## ----tcgaLoad-----------------------------------------------------------------
# By default MAF from MC3 project will be loaded
laml_mc3 = tcgaLoad(study = "LAML")
laml_mc3

# Change the source to Firehose
laml_fh = tcgaLoad(study = "LAML", source = "Firehose")
laml_fh

## ----eval=FALSE---------------------------------------------------------------
#  BiocManager::install(pkgs = "PoisonAlien/TCGAmutations")

## ----maf2mae------------------------------------------------------------------
laml_mae = maf2mae(m = laml)
laml_mae

## -----------------------------------------------------------------------------
sessionInfo()

