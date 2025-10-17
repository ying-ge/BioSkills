## ----echo = FALSE-------------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----code, echo = FALSE-------------------------------------------------------
code <- function(...) {
    cat(paste(..., sep = "\n"))
}

date = "`r doc_date()`"
pkg = "`r pkg_ver('BiocStyle')`"

## ----global_options, echo=FALSE-----------------------------------------------
short=TRUE #if short==TRUE, do not echo code chunks
debug=FALSE
knitr::opts_chunk$set(echo=!short, warning=debug, message=debug, error=FALSE,
               cache.path = "cache/",
               fig.path = "figures/")

## ----PFMatrix, echo=TRUE, eval=TRUE-------------------------------------------
library(TFBSTools)
## PFMatrix construction; Not all of the slots need to be initialised.
pfm <- PFMatrix(ID="MA0004.1", name="Arnt", 
                matrixClass="Zipper-Type", strand="+",
                bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                tags=list(family="Helix-Loop-Helix", species="10090",
                          tax_group="vertebrates",medline="7592839", 
                          type="SELEX",ACC="P53762", pazar_tf_id="TF0000003",
                          TFBSshape_ID="11", TFencyclopedia_ID="580"),
                profileMatrix=matrix(c(4L,  19L, 0L,  0L,  0L,  0L,
                                       16L, 0L,  20L, 0L,  0L,  0L,
                                       0L,  1L,  0L,  20L, 0L,  20L,
                                       0L,  0L,  0L,  0L,  20L, 0L),
                                     byrow=TRUE, nrow=4,
                                     dimnames=list(c("A", "C", "G", "T"))
                                     )
                )

pfm

## coerced to matrix
as.matrix(pfm)

## access the slots of pfm
ID(pfm)
name(pfm)
Matrix(pfm)
ncol(pfm)
length(pfm)

## convert a PFM to PWM, ICM
pwm <- toPWM(pfm, type="log2probratio", pseudocounts=0.8,
             bg=c(A=0.25, C=0.25, G=0.25, T=0.25))

icm <- toICM(pfm, pseudocounts=sqrt(rowSums(pfm)[1]), schneider=FALSE,
             bg=c(A=0.25, C=0.25, G=0.25, T=0.25))

## get the reverse complment matrix with all the same information except the strand.
pwmRevComp <- reverseComplement(pwm)

## ----PFMatrixList, echo=TRUE, eval=TRUE---------------------------------------
pfm2 <- pfm
pfmList <- PFMatrixList(pfm1=pfm, pfm2=pfm2, use.names=TRUE)
pfmList
names(pfmList)

## ----TFFMRead, echo=TRUE, eval=TRUE-------------------------------------------
  xmlFirst <- file.path(system.file("extdata", package="TFBSTools"),
                        "tffm_first_order.xml")
  tffmFirst <- readXMLTFFM(xmlFirst, type="First")
  tffm <- getPosProb(tffmFirst)

  xmlDetail <- file.path(system.file("extdata", package="TFBSTools"),
                         "tffm_detailed.xml")
  tffmDetail <- readXMLTFFM(xmlDetail, type="Detail")
  getPosProb(tffmDetail)

## ----searchDB, echo=TRUE, eval=TRUE-------------------------------------------
suppressMessages(library(JASPAR2014))
opts <- list()
opts[["species"]] <- 9606
opts[["name"]] <- "RUNX1"
opts[["type"]] <- "SELEX"
opts[["all_versions"]] <- TRUE
PFMatrixList <- getMatrixSet(JASPAR2014, opts)
PFMatrixList

opts2 <- list()
opts2[["type"]] <- "SELEX"
PFMatrixList2 <- getMatrixSet(JASPAR2014, opts2)
PFMatrixList2

## ----operateDb, echo=TRUE, eval=TRUE------------------------------------------
db <- "myMatrixDb.sqlite"
initializeJASPARDB(db, version="2014")
data("MA0043")
storeMatrix(db, MA0043)
deleteMatrixHavingID(db,"MA0043.1")
file.remove(db)

## ----PWMmatrixMethods, echo=TRUE, eval=TRUE-----------------------------------
pwm <- toPWM(pfm, pseudocounts=0.8)
pwm

## ----ICMmatrixMethods, echo=TRUE, eval=TRUE-----------------------------------
icm <- toICM(pfm, pseudocounts=0.8, schneider=TRUE)
icm

## ----seqLogo1, echo=TRUE, eval=TRUE, fig.width=6, fig.height=4----------------
seqLogo(icm)

## ----PFMSimi, echo=TRUE, eval=TRUE--------------------------------------------
## one to one comparison
data(MA0003.2)
data(MA0004.1)
pfmSubject <- MA0003.2
pfmQuery <-  MA0004.1
PFMSimilarity(pfmSubject, pfmQuery)

## one to several comparsion
PFMSimilarity(pfmList, pfmQuery)

## align IUPAC string
IUPACString <- "ACGTMRWSYKVHDBN"
PFMSimilarity(pfmList, IUPACString)

## ----PWMSimilarity, echo=TRUE, eval=TRUE--------------------------------------
data(MA0003.2)
data(MA0004.1)
pwm1 <- toPWM(MA0003.2, type="prob")
pwm2 <- toPWM(MA0004.1, type="prob")
PWMSimilarity(pwm1, pwm2, method="Euclidean")
PWMSimilarity(pwm1, pwm2, method="Pearson")
PWMSimilarity(pwm1, pwm2, method="KL")

## ----permuteMatrix, echo=TRUE, eval=TRUE--------------------------------------
## Matrice permutation
permuteMatrix(pfmQuery)
permuteMatrix(pfmList, type="intra")
permuteMatrix(pfmList, type="inter")

## ----samplingMatrix, echo=TRUE, eval=FALSE------------------------------------
#  ## Dirichlet model training
#  data(MA0003.2)
#  data(MA0004.1)
#  pfmList <- PFMatrixList(pfm1=MA0003.2, pfm2=MA0004.1, use.names=TRUE)
#  dmmParameters <- dmmEM(pfmList, K=6, alg="C")
#  ## Matrice sampling from trained Dirichlet model
#  pwmSampled <- rPWMDmm(MA0003.2, dmmParameters$alpha0, dmmParameters$pmix,
#                        N=1, W=6)

## ----TFFMFirstseqLogo, echo=TRUE, eval=TRUE, fig.width=6, fig.height=10-------
  ## sequence logo for First-order TFFM
  seqLogo(tffmFirst)

## ----TFFMDetailseqLogo, echo=TRUE, eval=TRUE, fig.width=6, fig.height=10------
  ## sequence logo for detailed TFFM
  seqLogo(tffmDetail)

## ----searchSeq, echo=TRUE, eval=TRUE------------------------------------------
library(Biostrings)
data(MA0003.2)
data(MA0004.1)
pwmList <- PWMatrixList(MA0003.2=toPWM(MA0003.2), MA0004.1=toPWM(MA0004.1),
                        use.names=TRUE)
subject <- DNAString("GAATTCTCTCTTGTTGTAGTCTCTTGACAAAATG")
siteset <- searchSeq(pwm, subject, seqname="seq1", min.score="60%", strand="*")

sitesetList <- searchSeq(pwmList, subject, seqname="seq1",
                         min.score="60%", strand="*")


## generate gff2 or gff3 style output
head(writeGFF3(siteset))
head(writeGFF3(sitesetList))
head(writeGFF2(siteset))

## get the relative scores
relScore(siteset)
relScore(sitesetList)

## calculate the empirical p-values of the scores
pvalues(siteset, type="TFMPvalue")
pvalues(siteset, type="sampling")

## ----searchAln, echo=TRUE, eval=TRUE------------------------------------------
library(Biostrings)
data(MA0003.2)
pwm <- toPWM(MA0003.2)
aln1 <- DNAString("ACTTCACCAGCTCCCTGGCGGTAAGTTGATC---AAAGG---AAACGCAAAGTTTTCAAG")
aln2 <- DNAString("GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC")
sitePairSet <- searchAln(pwm, aln1, aln2, seqname1="seq1", seqname2="seq2",
                         min.score="50%", cutoff=0.5,
                         strand="*", type="any")
## generate gff style output
head(writeGFF3(sitePairSet))
head(writeGFF2(sitePairSet))

## search the Axt alignment
library(CNEr)
axtFilesHg19DanRer7 <- file.path(system.file("extdata", package="TFBSTools"),
                                 "hg19.danRer7.net.axt")
axtHg19DanRer7 <- readAxt(axtFilesHg19DanRer7)
sitePairSet <-  searchAln(pwm, axtHg19DanRer7, min.score="80%",
                          windowSize=51L, cutoff=0.7, strand="*",
                          type="any", conservation=NULL, mc.cores=1)
GRangesTFBS <- toGRangesList(sitePairSet, axtHg19DanRer7)
GRangesTFBS$targetTFBS
GRangesTFBS$queryTFBS

## ----searchBSgenome, echo=TRUE, eval=FALSE------------------------------------
#  library(rtracklayer)
#  library(JASPAR2014)
#  library(BSgenome.Hsapiens.UCSC.hg19)
#  library(BSgenome.Mmusculus.UCSC.mm10)
#  pfm <- getMatrixByID(JASPAR2014, ID="MA0004.1")
#  pwm <- toPWM(pfm)
#  chain <- import.chain("Downloads/hg19ToMm10.over.chain")
#  sitePairSet <- searchPairBSgenome(pwm, BSgenome.Hsapiens.UCSC.hg19,
#                                   BSgenome.Mmusculus.UCSC.mm10,
#                                   chr1="chr1", chr2="chr1",
#                                   min.score="90%", strand="+", chain=chain)

## ----MEME-wrapper, echo=TRUE, eval=FALSE--------------------------------------
#  motifSet <- runMEME(file.path(system.file("extdata",
#                                            package="TFBSTools"), "crp0.s"),
#                      binary="meme",
#                      arguments=list("-nmotifs"=3)
#                     )
#  ## Get the sites sequences and surrounding sequences
#  sitesSeq(motifSet, type="all")
#  ## Get the sites sequences only
#  sitesSeq(motifSet, type="none")
#  consensusMatrix(motifSet)

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

