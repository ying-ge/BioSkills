## ----LoadPackageToDetermineVersion,echo=FALSE,message=FALSE,results='hide'----
options(width=65)
set.seed(0)
library(msa)
library(seqinr)
msaVersion <- packageDescription("msa")$Version
msaDateRaw <- packageDescription("msa")$Date
msaDateYear <- as.numeric(substr(msaDateRaw, 1, 4))
msaDateMonth <- as.numeric(substr(msaDateRaw, 6, 7))
msaDateDay <- as.numeric(substr(msaDateRaw, 9, 10))
msaDate <- paste(month.name[msaDateMonth], " ",
                 msaDateDay, ", ",
                 msaDateYear, sep="")

## ----InstallMSA,eval=FALSE-------------------------------------
#  if (!requireNamespace("BiocManager", quietly=TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("msa")

## ----LoadMSA,eval=FALSE----------------------------------------
#  library(msa)

## ----locateTeXshadeSty,eval=FALSE------------------------------
#  system.file("tex", "texshade.sty", package="msa")

## ----SimpleExFileNames-----------------------------------------
mySequenceFile <- system.file("examples", "exampleAA.fasta", package="msa")
mySequences <- readAAStringSet(mySequenceFile)
mySequences

## ----doAlignment-----------------------------------------------
myFirstAlignment <- msa(mySequences)
myFirstAlignment

## ----showWholeWidth--------------------------------------------
print(myFirstAlignment, show="complete")

## ----IntegratePDF2,eval=FALSE----------------------------------
#  msaPrettyPrint(myFirstAlignment, output="pdf", showNames="none",
#                 showLogo="none", askForOverwrite=FALSE, verbose=FALSE)

## ----VisualizePDF,results='asis'-------------------------------
msaPrettyPrint(myFirstAlignment, y=c(164, 213), output="asis",
               showNames="none", showLogo="none", askForOverwrite=FALSE)

## ----OtherAlgorithms-------------------------------------------
myClustalWAlignment <- msa(mySequences, "ClustalW")
myClustalWAlignment
myClustalOmegaAlignment <- msa(mySequences, "ClustalOmega")
myClustalOmegaAlignment
myMuscleAlignment <- msa(mySequences, "Muscle")
myMuscleAlignment

## ----helpPrint,eval=FALSE--------------------------------------
#  help("print,MsaDNAMultipleAlignment-method")

## ----printExamples---------------------------------------------
print(myFirstAlignment)
print(myFirstAlignment, show="complete")
print(myFirstAlignment, showConsensus=FALSE, halfNrow=3)
print(myFirstAlignment, showNames=FALSE, show="complete")

## ----maskExample-----------------------------------------------
myMaskedAlignment <- myFirstAlignment
colM <- IRanges(start=1, end=100)
colmask(myMaskedAlignment) <- colM
myMaskedAlignment

## ----unmaskedExample-------------------------------------------
unmasked(myMaskedAlignment)

## ----consensusExample1-----------------------------------------
conMat <- consensusMatrix(myFirstAlignment)
dim(conMat)
conMat[, 101:110]

## ----consensusExample2-----------------------------------------
conMat <- consensusMatrix(myMaskedAlignment)
conMat[, 95:104]

## ----consensusExample3-----------------------------------------
printSplitString <- function(x, width=getOption("width") - 1)
{
    starts <- seq(from=1, to=nchar(x), by=width)

    for (i in 1:length(starts))
        cat(substr(x, starts[i], starts[i] + width - 1), "\n")
}

printSplitString(msaConsensusSequence(myFirstAlignment))

## ----consensusExample4-----------------------------------------
printSplitString(msaConsensusSequence(myFirstAlignment, type="upperlower",
                                      thresh=c(40, 20)))

## ----consensusExample5-----------------------------------------
printSplitString(msaConsensusSequence(myMaskedAlignment, type="upperlower",
                                      thresh=c(40, 20)))

## ----conservationExample1--------------------------------------
data(BLOSUM62)
msaConservationScore(myFirstAlignment, BLOSUM62)

## ----conservationExample2--------------------------------------
msaConservationScore(myFirstAlignment, BLOSUM62, gapVsGap=0,
                     type="upperlower", thresh=c(40, 20))

## ----conservationExample3--------------------------------------
msaConservationScore(myMaskedAlignment, BLOSUM62, gapVsGap=0,
                     type="upperlower", thresh=c(40, 20))

## ----Hemoglobin1-----------------------------------------------
hemoSeq <- readAAStringSet(system.file("examples/HemoglobinAA.fasta",
                                       package="msa"))
hemoAln <- msa(hemoSeq)
hemoAln
hemoAln2 <- msaConvert(hemoAln, type="seqinr::alignment")

## ----Hemoglobin2-----------------------------------------------
library(seqinr)

d <- dist.alignment(hemoAln2, "identity")
as.matrix(d)[2:5, "HBA1_Homo_sapiens", drop=FALSE]

## ----HemoglobinTree,output.width='0.8\\textwidth',output.height='0.5\\textwidth',message=FALSE,results='hide'----
library(ape)

hemoTree <- nj(d)
plot(hemoTree, main="Phylogenetic Tree of Hemoglobin Alpha Sequences")

## ----Hemoglobin3-----------------------------------------------
hemoAln3 <- msaConvert(hemoAln, type="bios2mds::align")
str(hemoAln3)

## ----Hemoglobin4-----------------------------------------------
hemoAln4 <- as(hemoAln, "BStringSet")
hemoAln4

## ----ShowConsensusBottom,results="asis"------------------------
msaPrettyPrint(myFirstAlignment, output="asis", y=c(164, 213),
               subset=c(1:6), showNames="none", showLogo="none",
               consensusColor="ColdHot", showLegend=FALSE,
               askForOverwrite=FALSE)

## ----ShowLogoDefault,results="asis"----------------------------
msaPrettyPrint(myFirstAlignment, output="asis", y=c(164, 213),
               subset=c(1:6), showNames="none", showLogo="top",
               logoColors="rasmol", shadingMode="similar",
               showLegend=FALSE, askForOverwrite=FALSE)

## ----Shading1,results='asis'-----------------------------------
msaPrettyPrint(myFirstAlignment, output="asis", y=c(164, 213),
               showNames="none", shadingMode="similar",
               shadingColors="blues", showLogo="none",
               showLegend=FALSE, askForOverwrite=FALSE)

## ----Shading2,results='asis'-----------------------------------
msaPrettyPrint(myFirstAlignment, output="asis", y=c(164, 213),
               showNames="none", shadingMode="functional",
               shadingModeArg="structure",
               askForOverwrite=FALSE)

## ----ShowConsensusBottom2,results="asis"-----------------------
msaPrettyPrint(myFirstAlignment, output="asis", y=c(164, 213),
               subset=c(1:6), showNames="none", showLogo="none",
               consensusColor="ColdHot", showLegend=FALSE,
               shadingMode="similar", askForOverwrite=FALSE,
               furtherCode=c("\\defconsensus{.}{lower}{upper}",
                             "\\showruler{1}{top}"))

## ----SplitAlignmentIntoJunks,eval=FALSE------------------------
#  chunkSize <- 300 ## how much fits on one page depends on the length of
#                   ## names and the number of sequences;
#                   ## change to what suits your needs
#  
#  for (start in seq(1, ncol(aln), by=chunkSize))
#  {
#      end <- min(start + chunkSize - 1, ncol(aln))
#      alnPart <- DNAMultipleAlignment(subseq(unmasked(aln), start, end))
#  
#      msaPrettyPrint(x=alnPart, output="pdf", subset=NULL,
#                     file=paste0("aln_", start, "-", end, ".pdf"))
#  }

## ----GetBibTeX,eval=FALSE--------------------------------------
#  toBibtex(citation("msa"))

