## ----citation,eval=FALSE------------------------------------------------------
#  toBibtex(citation("minfi"))

## ----dependencies, warning=FALSE, message=FALSE-------------------------------
library(minfi)
library(minfiData)

## ----RGsetEx------------------------------------------------------------------
RGsetEx
## RGsetEx: RGChannelSet, 622,399 features
MsetEx <- preprocessRaw(RGsetEx)
## MsetEx: MethylSet, 485,512 features
GMsetEx <- mapToGenome(MsetEx)
## GMsetEx: GenomicMethylSet, 485,512 features

## ----baseDir------------------------------------------------------------------
baseDir <- system.file("extdata", package = "minfiData")
list.files(baseDir)

## ----baseDir2-----------------------------------------------------------------
list.files(file.path(baseDir, "5723646052"))

## ----sheet--------------------------------------------------------------------
targets <- read.metharray.sheet(baseDir)
targets

## ----BasenameColumn>----------------------------------------------------------
sub(baseDir, "", targets$Basename)

## ----readingTargets-----------------------------------------------------------
RGset <- read.metharray.exp(targets = targets)

## ----pData--------------------------------------------------------------------
RGset
pd <- pData(RGset)
pd[,1:4]

## ----read2--------------------------------------------------------------------
RGset2 <- read.metharray.exp(file.path(baseDir, "5723646052"))
RGset3 <- read.metharray.exp(baseDir, recursive = TRUE)

## ----sampleSheet2-------------------------------------------------------------
targets2 <- read.csv(file.path(baseDir, "SampleSheet.csv"), 
                     stringsAsFactors = FALSE, skip = 7)
targets2

## ----Basename-----------------------------------------------------------------
targets2$Basename <- file.path(baseDir, targets2$Sentrix_ID, 
                               paste0(targets2$Sentrix_ID, 
                                      targets2$Sentrix_Position))

## ----annotation---------------------------------------------------------------
annotation(RGsetEx)

## ----sessionInfo, results='asis', echo=FALSE----------------------------------
toLatex(sessionInfo())

