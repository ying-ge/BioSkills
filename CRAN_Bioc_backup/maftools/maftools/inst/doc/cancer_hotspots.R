## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(maftools)

## ----eval=FALSE---------------------------------------------------------------
#  can_hs_tbl = maftools::cancerhotspots(
#    bam = "Tumor.bam",
#    refbuild = "GRCh37",
#    mapq = 10,
#    sam_flag = 1024
#  )

## ----eval=FALSE---------------------------------------------------------------
#  head(can_hs_tbl)
#  
#  # loci fa_ref NT_change Hugo_Symbol Variant_Classification AA_change              Meta VAF A   T  G  C Ins Del
#  # 1: 1:2491289     NA       G>A    TNFRSF14      Missense_Mutation     C111Y    deleterious(0)   0 0   0 21  0   0   0
#  # 2: 1:2491290     NA       C>G    TNFRSF14      Missense_Mutation     C111W    deleterious(0)   0 0   0  0 21   0   0
#  # 3: 1:8073432     NA       T>G      ERRFI1      Missense_Mutation     K409N    deleterious(0)   0 1  64  0  1   0   0
#  # 4: 1:8073434     NA       T>G      ERRFI1      Missense_Mutation     K409Q deleterious(0.04)   0 0  63  0  0   0   0
#  # 5: 1:8074313     NA       T>A      ERRFI1      Nonsense_Mutation     K116*                     0 0 106  0  0   0   0
#  # 6: 1:9779982     NA       T>C      PIK3CD      Missense_Mutation     C416R   tolerated(0.26)   0 1  18  0  0   0   0

## -----------------------------------------------------------------------------
#Generate a sample loci - first two columns must contain chromosome name and position 
loci = data.table::data.table(chr = c("seq1", "seq2"), pos = c(1340, 1483))
loci

## ----eval=FALSE---------------------------------------------------------------
#  #Example BAM file from Rsamtools package
#  #By default position are assumed to be in 1-based coordinate system
#  bamfile = system.file("extdata", "ex1.bam", package = "Rsamtools")
#  loci_rc = maftools::bamreadcounts(bam = bamfile, loci = loci)
#  
#  loci_rc
#  # $ex1
#  #         loci fa_ref A  T G  C Ins Del
#  # 1: seq1:1340     NA 1  0 0 62   0   0
#  # 2: seq2:1483     NA 0 13 0  0   0   0

## -----------------------------------------------------------------------------
sessionInfo()

