## ----eval=TRUE----------------------------------------------------------------
library(DNAshapeR)

## ----eval=TRUE----------------------------------------------------------------
library(DNAshapeR)
fn <- system.file("extdata", "CGRsample.fa", package = "DNAshapeR")
pred <- getShape(fn)

## ----eval=FALSE---------------------------------------------------------------
#  # Install Bioconductor packages
#  if (!requireNamespace("BiocManager", quietly=TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("BSgenome.Scerevisiae.UCSC.sacCer3")
#  
#  library(BSgenome.Scerevisiae.UCSC.sacCer3)
#  
#  # Create a query GRanges object
#  gr <- GRanges(seqnames = c("chrI"),
#              strand = c("+", "-", "+"),
#              ranges = IRanges(start = c(100, 200, 300), width = 100))
#  getFasta(gr, Scerevisiae, width = 100, filename = "tmp.fa")
#  fn <- "tmp.fa"
#  pred <- getShape(fn)

## ----eval=FALSE---------------------------------------------------------------
#  # Install Bioconductor packages
#  library(BSgenome.Hsapiens.UCSC.hg19)
#  library(AnnotationHub)
#  
#  ah <- AnnotationHub()
#  ah <- subset(ah, species=="Homo sapiens")
#  ah <- query(ah, c("H3K4me3", "Gm12878", "Roadmap"))
#  getFasta(ah[[1]], Hsapiens, width = 150, filename = "tmp.fa")
#  fn <- "tmp.fa"
#  pred <- getShape(fn)

## ----eval=TRUE----------------------------------------------------------------
library(DNAshapeR)
fn_methy <- system.file("extdata", "MethylSample.fa", package = "DNAshapeR")
pred_methy <- getShape(fn_methy, methylate = TRUE)
pred_methy$MGW

## ----eval=TRUE----------------------------------------------------------------
library(DNAshapeR)
fn_methy <- system.file("extdata", "SingleSeqsample.fa", package = "DNAshapeR")
fn_methy_pos <- system.file("extdata", "MethylSamplePos.fa", package = "DNAshapeR")
pred_methy <- getShape(fn_methy, methylate = TRUE, methylatedPosFile = fn_methy_pos)
pred_methy$MGW

## ----fig.width=7, fig.height=7, fig.align='center', eval=TRUE-----------------
plotShape(pred$MGW)
#plotShape(pred$ProT)
#plotShape(pred$Roll)
#plotShape(pred$HelT)

## ----fig.width=7, fig.height=7, fig.align='center', eval=TRUE-----------------
library(fields)
heatShape(pred$ProT, 20)
#heatShape(pred$MGW, 20)
#heatShape(pred$Roll[1:500, 1:1980], 20)
#heatShape(pred$HelT[1:500, 1:1980], 20)

## ----fig.width=7, fig.height=7, fig.align='center', eval=TRUE-----------------
fn2 <- system.file("extdata", "SingleSeqsample.fa", package = "DNAshapeR")
pred2 <- getShape(fn2)
trackShape(fn2, pred2) # Only for single sequence file

## ----eval=TRUE----------------------------------------------------------------
library(Biostrings)
fn3 <- system.file("extdata", "PBMsample_short.fa", package = "DNAshapeR")
pred3 <- getShape(fn3)
featureType <- c("1-mer", "1-shape")
featureVector <- encodeSeqShape(fn3, pred3, featureType)
head(featureVector)

## ----eval=TRUE----------------------------------------------------------------
fn4 <- system.file("extdata", "PBMsample_short.s", package = "DNAshapeR")
experimentalData <- read.table(fn4)
df <- data.frame(affinity=experimentalData$V1, featureVector)

## ----eval=TRUE----------------------------------------------------------------
library(caret)

trainControl <- trainControl(method = "cv", number = 3, 
                savePredictions = TRUE)
model <- train (affinity~ ., data = df, 
                trControl=trainControl, method="lm", preProcess=NULL)
model

## ----eval=TRUE----------------------------------------------------------------
sessionInfo()

