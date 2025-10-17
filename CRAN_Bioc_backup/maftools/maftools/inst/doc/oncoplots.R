## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = 'left',
  fig.height = 5,
  fig.width = 10
)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(maftools)

## -----------------------------------------------------------------------------
#path to TCGA LAML MAF file
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')
#clinical information containing survival information and histology. This is optional
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools')

laml = read.maf(maf = laml.maf,
                clinicalData = laml.clin,
                verbose = FALSE)

## -----------------------------------------------------------------------------
#By default the function plots top20 mutated genes
oncoplot(maf = laml, draw_titv = TRUE)

## -----------------------------------------------------------------------------
#One can use any colors, here in this example color palette from RColorBrewer package is used
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)

print(vc_cols)

oncoplot(maf = laml, colors = vc_cols, top = 10)

## ----fig.height=5,fig.width=10, fig.align='left'------------------------------
#GISTIC results LAML
all.lesions =
  system.file("extdata", "all_lesions.conf_99.txt", package = "maftools")
amp.genes =
  system.file("extdata", "amp_genes.conf_99.txt", package = "maftools")
del.genes =
  system.file("extdata", "del_genes.conf_99.txt", package = "maftools")
scores.gis =
  system.file("extdata", "scores.gistic", package = "maftools")

#Read GISTIC results along with MAF
laml.plus.gistic = read.maf(
  maf = laml.maf,
  gisticAllLesionsFile = all.lesions,
  gisticAmpGenesFile = amp.genes,
  gisticDelGenesFile = del.genes,
  gisticScoresFile = scores.gis,
  isTCGA = TRUE,
  verbose = FALSE, 
  clinicalData = laml.clin
)

## ----fig.align='left',fig.height=5,fig.width=10, eval=T, fig.align='left'-----
oncoplot(maf = laml.plus.gistic, top = 10)

## -----------------------------------------------------------------------------
set.seed(seed = 1024)
barcodes = as.character(getSampleSummary(x = laml)[,Tumor_Sample_Barcode])
#Random 20 samples
dummy.samples = sample(x = barcodes,
                       size = 20,
                       replace = FALSE)

#Genarate random CN status for above samples
cn.status = sample(
  x = c('ShallowAmp', 'DeepDel', 'Del', 'Amp'),
  size = length(dummy.samples),
  replace = TRUE
)

custom.cn.data = data.frame(
  Gene = "DNMT3A",
  Sample_name = dummy.samples,
  CN = cn.status,
  stringsAsFactors = FALSE
)

head(custom.cn.data)

laml.plus.cn = read.maf(maf = laml.maf,
                        cnTable = custom.cn.data,
                        verbose = FALSE)

oncoplot(maf = laml.plus.cn, top = 5)

## ----fig.height=7,fig.width=10, eval=T, fig.align='left'----------------------
#Selected AML driver genes
aml_genes = c("TP53", "WT1", "PHF6", "DNMT3A", "DNMT3B", "TET1", "TET2", "IDH1", "IDH2", "FLT3", "KIT", "KRAS", "NRAS", "RUNX1", "CEBPA", "ASXL1", "EZH2", "KDM6A")

#Variant allele frequcnies (Right bar plot)
aml_genes_vaf = subsetMaf(maf = laml, genes = aml_genes, fields = "i_TumorVAF_WU", mafObj = FALSE)[,mean(i_TumorVAF_WU, na.rm = TRUE), Hugo_Symbol]
colnames(aml_genes_vaf)[2] = "VAF"
head(aml_genes_vaf)

#MutSig results (Right bar plot)
laml.mutsig = system.file("extdata", "LAML_sig_genes.txt.gz", package = "maftools")
laml.mutsig = data.table::fread(input = laml.mutsig)[,.(gene, q)]
laml.mutsig[,q := -log10(q)] #transoform to log10
head(laml.mutsig)


oncoplot(
  maf = laml,
  genes = aml_genes,
  leftBarData = aml_genes_vaf,
  leftBarLims = c(0, 100),
  rightBarData = laml.mutsig,
  rightBarLims = c(0, 20)
)

## -----------------------------------------------------------------------------
getClinicalData(x = laml)

## -----------------------------------------------------------------------------
oncoplot(maf = laml, genes = aml_genes, clinicalFeatures = 'FAB_classification')

## -----------------------------------------------------------------------------
#Color coding for FAB classification
fabcolors = RColorBrewer::brewer.pal(n = 8,name = 'Spectral')
names(fabcolors) = c("M0", "M1", "M2", "M3", "M4", "M5", "M6", "M7")
fabcolors = list(FAB_classification = fabcolors)

print(fabcolors)

oncoplot(
  maf = laml, genes = aml_genes,
  clinicalFeatures = 'FAB_classification',
  sortByAnnotation = TRUE,
  annotationColor = fabcolors
)

## -----------------------------------------------------------------------------
oncoplot(maf = laml, genes = aml_genes,
         additionalFeature = c("Tumor_Seq_Allele2", "C"))

## -----------------------------------------------------------------------------
getFields(x = laml)

## ----fig.height = 8, fig.width = 10-------------------------------------------
oncoplot(maf = laml, pathways = "sigpw", gene_mar = 8, fontSize = 0.6, topPathways = 5)

## ----fig.height = 8, fig.width = 10-------------------------------------------
oncoplot(maf = laml, pathways = "smgbp", gene_mar = 8, fontSize = 0.8, topPathways = 5)

## ----fig.height = 7, fig.width = 10-------------------------------------------
pathways = data.frame(
  Genes = c(
    "TP53",
    "WT1",
    "PHF6",
    "DNMT3A",
    "DNMT3B",
    "TET1",
    "TET2",
    "IDH1",
    "IDH2",
    "FLT3",
    "KIT",
    "KRAS",
    "NRAS",
    "RUNX1",
    "CEBPA",
    "ASXL1",
    "EZH2",
    "KDM6A"
  ),
  Pathway = rep(c(
    "TSG", "DNAm", "Signalling", "TFs", "ChromMod"
  ), c(3, 6, 4, 2, 3)),
  stringsAsFactors = FALSE
)

head(pathways)

oncoplot(maf = laml, pathways = pathways, gene_mar = 8, fontSize = 0.6)

## -----------------------------------------------------------------------------
oncoplot(maf = laml, pathways = "sigpw", gene_mar = 8, fontSize = 0.6, topPathways = 5, collapsePathway = TRUE)

## ----fig.height = 8, fig.width = 10-------------------------------------------
oncoplot(
  maf = laml.plus.gistic,
  draw_titv = TRUE,
  pathways = pathways,
  clinicalFeatures = c('FAB_classification', 'Overall_Survival_Status'),
  sortByAnnotation = TRUE,
  additionalFeature = c("Tumor_Seq_Allele2", "C"),
  leftBarData = aml_genes_vaf,
  leftBarLims = c(0, 100),
  rightBarData = laml.mutsig[,.(gene, q)],
)

## -----------------------------------------------------------------------------
sessionInfo()

