## ----message = FALSE----------------------------------------------------------
library(SICtools)

## -----------------------------------------------------------------------------
bam1 <- system.file(package='SICtools','extdata','example1.bam')
bam2 <- system.file(package='SICtools','extdata','example2.bam')
refFsa <- system.file(package='SICtools','extdata','example.ref.fasta')

snpDiffDf <- snpDiff(bam1,bam2,refFsa,'chr04',962501,1026983,pValueCutOff=1,baseDistCutOff=0)
snpDiffDf

## -----------------------------------------------------------------------------
snpDiffDfSort <- snpDiffDf[order(snpDiffDf$p.value,snpDiffDf$d.value),]
snpDiffDfSort

## -----------------------------------------------------------------------------
indelDiffDf <- indelDiff(bam1,bam2,refFsa,'chr07',828514,828914,pValueCutOff=1,gtDistCutOff=0)
indelDiffDfSort <- indelDiffDf[order(indelDiffDf$p.value,indelDiffDf$d.value),]
indelDiffDfSort

## -----------------------------------------------------------------------------
sessionInfo()

