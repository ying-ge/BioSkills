## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
# ## Install MSstatsTMT package from Bioconductor
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("MSstatsTMT")

library(MSstatsTMT)

## ----eval=T, echo=F, warning=F------------------------------------------------
head(input.pd)

## -----------------------------------------------------------------------------
# read in PD PSM sheet
# raw.pd <- read.delim("161117_SILAC_HeLa_UPS1_TMT10_5Mixtures_3TechRep_UPSdB_Multiconsensus_PD22_Intensity_PSMs.txt")
head(raw.pd)

# Read in annotation including condition and biological replicates per run and channel.
# Users should make this annotation file. It is not the output from Proteome Discoverer.
# annotation.pd <- read.csv(file="PD_Annotation.csv", header=TRUE)
head(annotation.pd)

# use Protein.Accessions as protein name
input.pd <- PDtoMSstatsTMTFormat(raw.pd, annotation.pd, 
                                 which.proteinid = "Protein.Accessions")
head(input.pd)

# use Master.Protein.Accessions as protein name
input.pd.master <- PDtoMSstatsTMTFormat(raw.pd, annotation.pd,
                                 which.proteinid = "Master.Protein.Accessions")
head(input.pd.master)

## -----------------------------------------------------------------------------
# Read in MaxQuant files
# proteinGroups <- read.table("proteinGroups.txt", sep="\t", header=TRUE)

# evidence <- read.table("evidence.txt", sep="\t", header=TRUE)

# Users should make this annotation file. It is not the output from MaxQuant.
# annotation.mq <- read.csv(file="MQ_Annotation.csv", header=TRUE)

input.mq <- MaxQtoMSstatsTMTFormat(evidence, proteinGroups, annotation.mq)
head(input.mq)

## -----------------------------------------------------------------------------
# Read in SpectroMine PSM report
# raw.mine <- read.csv('20180831_095547_CID-OT-MS3-Short_PSM Report_20180831_103118.xls', sep="\t")

# Users should make this annotation file. It is not the output from SpectroMine
# annotation.mine <- read.csv(file="Mine_Annotation.csv", header=TRUE)

input.mine <- SpectroMinetoMSstatsTMTFormat(raw.mine, annotation.mine)
head(input.mine)

## -----------------------------------------------------------------------------
# read in MSstatsTMT report from OpenMS
# raw.om <- read.csv("OpenMS_20200222/20200225_MSstatsTMT_OpenMS_Export.csv")
head(raw.om)

# the function only requries one input file
input.om <- OpenMStoMSstatsTMTFormat(raw.om)
head(input.om)

## -----------------------------------------------------------------------------
# Example code is skipped for Philosopher Converter 
# since the input is a path to the folder with all the Philosopher msstats csv files

## ----message=F, warning=F, results='hide'-------------------------------------
# use MSstats for protein summarization
quant.msstats <- proteinSummarization(input.pd,
                                      method="msstats",
                                      global_norm=TRUE,
                                      reference_norm=TRUE,
                                      remove_norm_channel = TRUE,
                                      remove_empty_channel = TRUE)

## -----------------------------------------------------------------------------
head(quant.pd.msstats$ProteinLevelData)

## ----message=F, warning=F, results='hide'-------------------------------------
# use Median for protein summarization
quant.median <- proteinSummarization(input.pd,
                                     method="Median",
                                     global_norm=TRUE,
                                     reference_norm=TRUE,
                                     remove_norm_channel = TRUE,
                                     remove_empty_channel = TRUE)

## -----------------------------------------------------------------------------
head(quant.median$ProteinLevelData)

## ----message=F, warning=F, results='hide'-------------------------------------
## Profile plot without norm channnels and empty channels
dataProcessPlotsTMT(data=quant.msstats,
                     type = 'ProfilePlot',
                     width = 21, # adjust the figure width since there are 15 TMT runs.
                     height = 7)

## ----message=F, warning=F, results='hide'-------------------------------------
dataProcessPlotsTMT(data=quant.msstats,
                    type='ProfilePlot', # choice of visualization
                    width = 21,
                    height = 7,
                    which.Protein = 'P04406') 

## -----------------------------------------------------------------------------
## Quality control plot 
# dataProcessPlotsTMT(data=quant.msstats, 
                     # type='QCPlot',
                     # width = 21, # adjust the figure width since there are 15 TMT runs. 
                     # height = 7)

## ----message=F, warning=F, results='hide'-------------------------------------
# test for all the possible pairs of conditions
test.pairwise <- groupComparisonTMT(quant.msstats, moderated = TRUE)

## -----------------------------------------------------------------------------
# Show test result
# Label : which comparison is used
# log2FC : estimated log2 fold change between two conditions (the contrast)
# adj.pvalue : adjusted p value
head(test.pairwise$ComparisonResult)

## -----------------------------------------------------------------------------
# Check the conditions in the protein level data
levels(quant.msstats$ProteinLevelData$Condition)
# Only compare condition 0.125 and 1
comparison<-matrix(c(-1,0,0,1),nrow=1)
# Set the names of each row
row.names(comparison)<-"1-0.125"
# Set the column names
colnames(comparison)<- c("0.125", "0.5", "0.667", "1")
comparison

## ----message=F, warning=F, results='hide'-------------------------------------
test.contrast <- groupComparisonTMT(data = quant.msstats, contrast.matrix = comparison, moderated = TRUE)

## -----------------------------------------------------------------------------
head(test.contrast$ComparisonResult)

