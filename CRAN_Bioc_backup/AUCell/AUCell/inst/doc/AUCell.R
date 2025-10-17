## ----echo=FALSE, message=FALSE, warning=FALSE---------------------------------
suppressPackageStartupMessages({
library(AUCell)
#library(Biobase)
library(GSEABase)
library(data.table)
library(DT)
library(NMF)
library(plotly)
library(GEOquery)
library(Matrix)
# library(doMC);library(doRNG) # Loaded by AUCell, to avoid messages. Not available in windows?
})

# To build a personalized report, update this working directory:
dir.create("AUCell_tutorial")
knitr::opts_knit$set(root.dir = 'AUCell_tutorial')

## ----Overview, eval=FALSE-----------------------------------------------------
#  library(AUCell)
#  geneSets <- list(geneSet1=c("gene1", "gene2", "gene3"))
#  
#  # Calculate enrichment scores
#  cells_AUC <- AUCell_run(exprMatrix, geneSets, aucMaxRank=nrow(cells_rankings)*0.05)
#  
#  # Optional: Set the assignment thresholds
#  par(mfrow=c(3,3))
#  set.seed(123)
#  cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)

## ----citation, echo=FALSE-----------------------------------------------------
print(citation("AUCell")[1], style="textVersion")

## ----setup, eval=FALSE--------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly=TRUE))
#      install.packages("BiocManager")
#  # To support paralell execution:
#  BiocManager::install(c("doMC", "doRNG","doSNOW"))
#  # For the main example:
#  BiocManager::install(c("mixtools", "SummarizedExperiment"))
#  # For the examples in the follow-up section of the tutorial:
#  BiocManager::install(c("DT", "plotly", "NMF", "d3heatmap", "shiny", "rbokeh",
#                         "dynamicTreeCut","R2HTML","Rtsne", "zoo"))

## ----vignette, eval=FALSE-----------------------------------------------------
#  # Explore tutorials in the web browser:
#  browseVignettes(package="AUCell")
#  
#  # Commnad line-based:
#  vignette(package="AUCell") # list
#  vignette("AUCell") # open

## ----editRmd, eval=FALSE------------------------------------------------------
#  vignetteFile <- paste(file.path(system.file('doc', package='AUCell')), "AUCell.Rmd", sep="/")
#  # Copy to edit as markdown
#  file.copy(vignetteFile, ".")
#  # Alternative: extract R code
#  Stangle(vignetteFile)

## ----setwd, eval=FALSE--------------------------------------------------------
#  dir.create("AUCell_tutorial")
#  setwd("AUCell_tutorial") # or in the first code chunk (kntr options), if running as Notebook

## ----loadingExprMat, eval=FALSE-----------------------------------------------
#  # i.e. Reading from a text file
#  exprMatrix <- read.table("myCountsMatrix.tsv")
#  
#  # or single-cell experiment
#  exprMatrix <- assay(mySingleCellExperiment)
#  
#  ### Convert to sparse:
#  exprMatrix <- as(exprMatrix, "dgCMatrix")

## ----GEOdataset, cache=TRUE, results='hide', message=FALSE, eval=TRUE---------
# (Downloads the data)
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
library(GEOquery)
# gse <- getGEO('GSE60361') # does not work, the matrix is in a suppl file
geoFile <- "GSE60361_C1-3005-Expression.txt.gz"
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60361/suppl/GSE60361_C1-3005-Expression.txt.gz", destfile = geoFile)

library(data.table)
exprMatrix <- fread(geoFile, sep="\t")
geneNames <- unname(unlist(exprMatrix[,1, with=FALSE]))
exprMatrix <- as.matrix(exprMatrix[,-1, with=FALSE])
rownames(exprMatrix) <- geneNames
exprMatrix <- exprMatrix[unique(rownames(exprMatrix)),]
dim(exprMatrix)
exprMatrix[1:5,1:4]

# Remove file(s) downloaded:
file.remove(geoFile)

# Convert to sparse
library(Matrix)
exprMatrix <- as(exprMatrix, "dgCMatrix")

# Save for future use
mouseBrainExprMatrix <- exprMatrix
save(mouseBrainExprMatrix, file="exprMatrix_MouseBrain.RData")

## ----randomSamples------------------------------------------------------------
# load("exprMatrix_MouseBrain.RData")
set.seed(333)
exprMatrix <- mouseBrainExprMatrix[sample(rownames(mouseBrainExprMatrix), 5000),]

## ----dimExprMat---------------------------------------------------------------
dim(exprMatrix)

## ----geneSetsFake-------------------------------------------------------------
library(GSEABase)
genes <- c("gene1", "gene2", "gene3")
geneSets <- GeneSet(genes, setName="geneSet1")
geneSets

## ----geneSets-----------------------------------------------------------------
library(AUCell)
library(GSEABase)
gmtFile <- paste(file.path(system.file('examples', package='AUCell')), "geneSignatures.gmt", sep="/")
geneSets <- getGmt(gmtFile)

## ----geneSetsNgenes-----------------------------------------------------------
geneSets <- subsetGeneSets(geneSets, rownames(exprMatrix)) 
cbind(nGenes(geneSets))

## ----geneSetsRename-----------------------------------------------------------
geneSets <- setGeneSetNames(geneSets, newNames=paste(names(geneSets), " (", nGenes(geneSets) ,"g)", sep=""))

## ----hkGs---------------------------------------------------------------------
# Random
set.seed(321)
extraGeneSets <- c(
  GeneSet(sample(rownames(exprMatrix), 50), setName="Random (50g)"),
  GeneSet(sample(rownames(exprMatrix), 500), setName="Random (500g)"))

countsPerGene <- apply(exprMatrix, 1, function(x) sum(x>0))
# Housekeeping-like
extraGeneSets <- c(extraGeneSets,
                   GeneSet(sample(names(countsPerGene)[which(countsPerGene>quantile(countsPerGene, probs=.95))], 100), setName="HK-like (100g)"))

geneSets <- GeneSetCollection(c(geneSets,extraGeneSets))
names(geneSets)

## ----runAUCell, cache=TRUE----------------------------------------------------
cells_AUC <- AUCell_run(exprMatrix, geneSets)
save(cells_AUC, file="cells_AUC.RData")

## -----------------------------------------------------------------------------
cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=FALSE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)

## ----buildRankings, cache=TRUE, fig.width=5, fig.height=5---------------------
cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=TRUE)

## -----------------------------------------------------------------------------
cells_rankings

## ----saveRankings, eval=FALSE-------------------------------------------------
#  save(cells_rankings, file="cells_rankings.RData")

## ----explainAUC, echo=FALSE---------------------------------------------------
geneSet <- geneSets[[1]]
gSetRanks <- cells_rankings[geneIds(geneSet),]

par(mfrow=c(1,2))
set.seed(222)
aucMaxRank <- nrow(cells_rankings)*0.05
na <- sapply(c(1, 2000), function(i){
  x <- sort(getRanking(gSetRanks[,i]))
  aucCurve <- cbind(c(0, x, nrow(cells_rankings)), c(0:length(x), length(x)))
  op <- par(mar=c(5, 6, 4, 2) + 0.1)
  plot(aucCurve, 
       type="s", col="darkblue", lwd=1, 
       xlab="Gene rank", ylab=paste("# genes in the gene set \n Gene set:", setName(geneSet)), 
       xlim=c(0, aucMaxRank*2), ylim=c(0, nGenes(geneSet)*.20), 
       main="Recovery curve", 
       sub=paste("Cell:", colnames(gSetRanks)[i]))
  aucShade <- aucCurve[which(aucCurve[,1] < aucMaxRank),]
  aucShade <- rbind(aucShade, c(aucMaxRank, nrow(aucShade)))
  aucShade[,1] <-  aucShade[,1]-1
  aucShade <- rbind(aucShade, c(max(aucShade),0))
  polygon(aucShade, col="#0066aa40", border=FALSE)
  
  abline(v=aucMaxRank, lty=2)
  text(aucMaxRank-50, 5, "AUC")
})

## ----calcAUC, cache=TRUE, warning=FALSE---------------------------------------
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
save(cells_AUC, file="cells_AUC.RData")

## ----exploreThresholds, warning=FALSE, fig.width=7, fig.height=7--------------
set.seed(333)
par(mfrow=c(3,3)) 
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 
warningMsg <- sapply(cells_assignment, function(x) x$aucThr$comment)
warningMsg[which(warningMsg!="")]

## ----explThr1-----------------------------------------------------------------
cells_assignment$Oligodendrocyte_Cahoy$aucThr$thresholds

## ----explThr2-----------------------------------------------------------------
cells_assignment$Oligodendrocyte_Cahoy$aucThr$selected
# getThresholdSelected(cells_assignment)

## ----cellsAssigned------------------------------------------------------------
oligodencrocytesAssigned <- cells_assignment$Oligodendrocyte_Cahoy$assignment
length(oligodencrocytesAssigned)
head(oligodencrocytesAssigned)

## ----AUCell_plot--------------------------------------------------------------
geneSetName <- rownames(cells_AUC)[grep("Oligodendrocyte_Cahoy", rownames(cells_AUC))]
AUCell_plotHist(cells_AUC[geneSetName,], aucThr=0.25)
abline(v=0.25)

## ----explThr3-----------------------------------------------------------------
newSelectedCells <- names(which(getAUC(cells_AUC)[geneSetName,]>0.08))
length(newSelectedCells)
head(newSelectedCells)

## ----explAssignment-----------------------------------------------------------
cellsAssigned <- lapply(cells_assignment, function(x) x$assignment)
assignmentTable <- reshape2::melt(cellsAssigned, value.name="cell")
colnames(assignmentTable)[2] <- "geneSet"
head(assignmentTable)

## ----assignmentMat------------------------------------------------------------
assignmentMat <- table(assignmentTable[,"geneSet"], assignmentTable[,"cell"])
assignmentMat[,1:2]

## ----assignHeatmap------------------------------------------------------------
set.seed(123)
miniAssigMat <- assignmentMat[,sample(1:ncol(assignmentMat),100)]
library(NMF)
aheatmap(miniAssigMat, scale="none", color="black", legend=FALSE)

## ----assignHeatmap_interactive, eval=FALSE------------------------------------
#  # if (!requireNamespace("d3heatmap", quietly = TRUE)) BiocManager::install("d3heatmap")
#  # d3heatmap(matrix(miniAssigMat), scale="none", colors=c("white", "black"))
#  
#  library(DT)
#  datatable(assignmentTable, options = list(pageLength = 10), filter="top")

## ----loadtSNE, fig.width=4, fig.height=4--------------------------------------
# Load the tSNE (included in the package)
load(paste(file.path(system.file('examples', package='AUCell')), "cellsTsne.RData", sep="/"))
cellsTsne <- cellsTsne$Y
plot(cellsTsne, pch=16, cex=.3)

## ----runTsne, eval=FALSE------------------------------------------------------
#  load("exprMatrix_AUCellVignette_MouseBrain.RData")
#  sumByGene <- apply(mouseBrainExprMatrix, 1, sum)
#  exprMatSubset <- mouseBrainExprMatrix[which(sumByGene>0),]
#  logMatrix <- log2(exprMatSubset+1)
#  logMatrix <- as.matrix(logMatrix)
#  
#  library(Rtsne)
#  set.seed(123)
#  cellsTsne <- Rtsne(t(logMatrix))
#  rownames(cellsTsne$Y) <- colnames(logMatrix)
#  colnames(cellsTsne$Y) <- c("tsne1", "tsne2")
#  save(cellsTsne, file="cellsTsne.RData")

## ----plotTsneCode, fig.width=7, fig.height=6----------------------------------
selectedThresholds <- getThresholdSelected(cells_assignment)
par(mfrow=c(2,3)) # Splits the plot into two rows and three columns
for(geneSetName in names(selectedThresholds)[1:6])
{
  nBreaks <- 5 # Number of levels in the color palettes
  # Color palette for the cells that do not pass the threshold
  colorPal_Neg <- grDevices::colorRampPalette(c("black","blue", "skyblue"))(nBreaks)
  # Color palette for the cells that pass the threshold
  colorPal_Pos <- grDevices::colorRampPalette(c("pink", "magenta", "red"))(nBreaks)
  
  # Split cells according to their AUC value for the gene set
  passThreshold <- getAUC(cells_AUC)[geneSetName,] >  selectedThresholds[geneSetName]
  if(sum(passThreshold) >0 )
  {
    aucSplit <- split(getAUC(cells_AUC)[geneSetName,], passThreshold)
    
    # Assign cell color
    cellColor <- c(setNames(colorPal_Neg[cut(aucSplit[[1]], breaks=nBreaks)], names(aucSplit[[1]])), 
    setNames(colorPal_Pos[cut(aucSplit[[2]], breaks=nBreaks)], names(aucSplit[[2]])))
    
    # Plot
    plot(cellsTsne, main=geneSetName,
    sub="Pink/red cells pass the threshold",
    col=cellColor[rownames(cellsTsne)], pch=16) 
  }
}

## ----tsneThreshold------------------------------------------------------------
selectedThresholds[2] <-  0.25
par(mfrow=c(2,3))
AUCell_plotTSNE(tSNE=cellsTsne, exprMat=exprMatrix, 
cellsAUC=cells_AUC[1:2,], thresholds=selectedThresholds)

## ----eval=FALSE---------------------------------------------------------------
#  library(shiny); library(rbokeh)
#  exprMat <- log2(mouseBrainExprMatrix+1) # (back to the full matrix)
#  
#  # Create app
#  aucellApp <- AUCell_createViewerApp(auc=cells_AUC, thresholds=selectedThresholds,
#  tSNE=cellsTsne, exprMat=exprMat)
#  
#  # Run (the exact commands depend on the R settings, see Shiny's doc for help)
#  options(shiny.host="0.0.0.0")
#  savedSelections <- runApp(aucellApp)
#  
#  # Other common settings:
#  # host = getOption("shiny.host", "127.0.0.1")
#  # shinyApp(ui=aucellApp$ui, server=aucellApp$server)
#  
#  # Optional: Visualize extra information about the cells (categorical variable)
#  cellInfoFile <- paste(file.path(system.file('examples', package='AUCell')),
#  "mouseBrain_cellLabels.tsv", sep="/")
#  cellInfo <- read.table(cellInfoFile, row.names=1, header=TRUE, sep="\t")
#  head(cellInfo)
#  
#  # To modify the cell assignment based on those thresholds:
#  newAssignments <- AUCell_assignCells(cells_AUC, thresholds=savedSelections$thresholds)
#  newAssignments <- getAssignments(newAssignments)
#  lengths(newAssignments)# number of cells assigned

## ----tSNE_interactive, eval=FALSE---------------------------------------------
#  geneSetName <- "Astrocyte_Cahoy (526g)"
#  
#  library(rbokeh)
#  tSNE.df <- data.frame(cellsTsne, cell=rownames(cellsTsne))
#  figure() %>%
#    ly_points(tsne1, tsne2, data=tSNE.df, size=1, legend = FALSE,
#      hover=cell, color=getAUC(cells_AUC)[geneSetName,rownames(tSNE.df)]) %>%
#    set_palette(continuous_color = pal_gradient(c("lightgrey", "pink", "red")))

## -----------------------------------------------------------------------------
logMat <- log2(exprMatrix+2)
meanByGs <- t(sapply(geneSets, function(gs) colMeans(logMat[geneIds(gs),])))
rownames(meanByGs) <- names(geneSets)

## ----meanTsne, fig.width=7, fig.height=8--------------------------------------
colorPal <- grDevices::colorRampPalette(c("black", "red"))(nBreaks)
par(mfrow=c(3,3))
for(geneSetName in names(geneSets))
{
  cellColor <- setNames(colorPal[cut(meanByGs[geneSetName,], breaks=nBreaks)], names(meanByGs[geneSetName,]))
  plot(cellsTsne, main=geneSetName, axes=FALSE, xlab="", ylab="",
  sub="Expression mean",
  col=cellColor[rownames(cellsTsne)], pch=16)
}

AUCell_plotTSNE(tSNE=cellsTsne, exprMat=exprMatrix, plots = "AUC",
cellsAUC=cells_AUC[1:9,], thresholds=selectedThresholds)

## ----fig.height=4, fig.width=3.5----------------------------------------------
nGenesPerCell <- apply(exprMatrix, 2, function(x) sum(x>0))
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))
cellColorNgenes <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(nGenesPerCell,breaks=10, right=FALSE,include.lowest=TRUE))], names(nGenesPerCell))
plot(cellsTsne, axes=FALSE, xlab="", ylab="",
sub="Number of detected genes",
col=cellColorNgenes[rownames(cellsTsne)], pch=16)

## -----------------------------------------------------------------------------
# "Real" cell type (e.g. provided in the publication)
cellLabels <- paste(file.path(system.file('examples', package='AUCell')), "mouseBrain_cellLabels.tsv", sep="/")
cellLabels <- read.table(cellLabels, row.names=1, header=TRUE, sep="\t")

## -----------------------------------------------------------------------------
# Confusion matrix:
cellTypeNames <- unique(cellLabels[,"level1class"])
confMatrix <- t(sapply(cells_assignment[c(1,4,2,5,3,6:9)], 
function(x) table(cellLabels[x$assignment,])[cellTypeNames]))
colnames(confMatrix) <- cellTypeNames
confMatrix[which(is.na(confMatrix), arr.ind=TRUE)] <- 0
confMatrix

## -----------------------------------------------------------------------------
date()
sessionInfo()

