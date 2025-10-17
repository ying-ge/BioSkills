## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(maftools)

## ----eval=FALSE---------------------------------------------------------------
#  remotes::install_github(repo = 'VanLoo-lab/ascat/ASCAT')

## ----eval=FALSE---------------------------------------------------------------
#  #Matched normal BAM files are strongly recommended
#  counts = maftools::gtMarkers(t_bam = "tumor.bam",
#                               n_bam = "normal.bam",
#                               build = "hg19")

## ----eval=FALSE---------------------------------------------------------------
#  library(ASCAT)
#  ascat.bc = maftools::prepAscat(t_counts = "tumor_nucleotide_counts.tsv",
#                                 n_counts = "normal_nucleotide_counts.tsv",
#                                 sample_name = "tumor")
#  
#  # Library sizes:
#  # Tumor:  1830168947
#  # Normal: 1321201848
#  # Library size difference: 1.385
#  # ------
#  # Counts file: tumor_nucleotide_counts.tsv
#  # Markers: 932148
#  # Removed 2982 duplicated loci
#  # Markers > 15: 928607
#  # ------
#  # Counts file: normal_nucleotide_counts.tsv
#  # Markers: 932148
#  # Removed 2982 duplicated loci
#  # Markers > 15: 928311
#  # ------
#  # Final number SNPs: 928107
#  # Generated following files:
#  # tumor_nucleotide_counts.tumour.BAF.txt
#  # tumor_nucleotide_counts.tumour.logR.txt
#  # tumor_nucleotide_counts.normal.BAF.txt
#  # tumor_nucleotide_counts.normal.logR.txt
#  # ------

## ----eval=FALSE---------------------------------------------------------------
#  
#  ascat.bc = ASCAT::ascat.loadData(
#    Tumor_LogR_file = "tumor_nucleotide_counts.tumour.logR.txt",
#    Tumor_BAF_file = "tumor_nucleotide_counts.tumour.BAF.txt",
#    Germline_LogR_file = "tumor_nucleotide_counts.normal.logR.txt",
#    Germline_BAF_file = "tumor_nucleotide_counts.normal.BAF.txt",
#    chrs = c(1:22, "X", "Y"),
#    sexchromosomes = c("X", "Y")
#  )
#  
#  ASCAT::ascat.plotRawData(ASCATobj = ascat.bc, img.prefix = "tumor")
#  ascat.bc = ASCAT::ascat.aspcf(ascat.bc)
#  ASCAT::ascat.plotSegmentedData(ascat.bc)
#  ascat.output = ASCAT::ascat.runAscat(ascat.bc)

## ----eval=FALSE---------------------------------------------------------------
#  ascat.bc = maftools::prepAscat_t(t_counts = "tumor_nucleotide_counts.tsv", sample_name = "tumor_only")
#  
#  # Library sizes:
#  # Tumor: 1830168947
#  # Counts file: tumor_nucleotide_counts.tsv
#  # Markers: 932148
#  # Removed 2982 duplicated loci
#  # Markers > 15: 928607
#  # Median depth of coverage (autosomes): 76
#  # ------
#  # Generated following files:
#  # tumor_only.tumour.BAF.txt
#  # tumor_only.tumour.logR.txt
#  # ------

## ----eval=FALSE---------------------------------------------------------------
#  ascat.bc = ASCAT::ascat.loadData(
#    Tumor_LogR_file = "tumor_only.tumour.logR.txt",
#    Tumor_BAF_file = "tumor_only.tumour.BAF.txt",
#    chrs = c(1:22, "X", "Y"),
#    sexchromosomes = c("X", "Y")
#  )
#  
#  ASCAT::ascat.plotRawData(ASCATobj = ascat.bc, img.prefix = "tumor_only")
#  ascat.gg = ASCAT::ascat.predictGermlineGenotypes(ascat.bc)
#  ascat.bc = ASCAT::ascat.aspcf(ascat.bc, ascat.gg=ascat.gg)
#  ASCAT::ascat.plotSegmentedData(ascat.bc)
#  ascat.output = ASCAT::ascat.runAscat(ascat.bc)

## ----eval=FALSE---------------------------------------------------------------
#  maftools::segmentLogR(tumor_logR = "tumor.tumour.logR.txt", sample_name = "tumor")
#  
#  # Analyzing: tumor
#  #   current chromosome: 1
#  #   current chromosome: 2
#  #   current chromosome: 3
#  #   current chromosome: 4
#  #   current chromosome: 5
#  #   current chromosome: 6
#  #   current chromosome: 7
#  #   current chromosome: 8
#  #   current chromosome: 9
#  #   current chromosome: 10
#  #   current chromosome: 11
#  #   current chromosome: 12
#  #   current chromosome: 13
#  #   current chromosome: 14
#  #   current chromosome: 15
#  #   current chromosome: 16
#  #   current chromosome: 17
#  #   current chromosome: 18
#  #   current chromosome: 19
#  #   current chromosome: 20
#  #   current chromosome: 21
#  #   current chromosome: 22
#  #   current chromosome: MT
#  #   current chromosome: X
#  #   current chromosome: Y
#  # Segments are written to: tumor_only.tumour_cbs.seg
#  # Segments are plotted to: tumor_only.tumour_cbs.png

## ----eval=FALSE---------------------------------------------------------------
#  plotMosdepth(
#    t_bed = "tumor.regions.bed.gz",
#    n_bed = "normal.regions.bed.gz",
#    segment = TRUE,
#    sample_name = "tumor"
#  )
#  
#  # Coverage ratio T/N: 1.821
#  # Running CBS segmentation:
#  # Analyzing: tumor01
#  #   current chromosome: 1
#  #   current chromosome: 2
#  #   current chromosome: 3
#  #   current chromosome: 4
#  #   current chromosome: 5
#  #   current chromosome: 6
#  #   current chromosome: 7
#  #   current chromosome: 8
#  #   current chromosome: 9
#  #   current chromosome: 10
#  #   current chromosome: 11
#  #   current chromosome: 12
#  #   current chromosome: 13
#  #   current chromosome: 14
#  #   current chromosome: 15
#  #   current chromosome: 16
#  #   current chromosome: 17
#  #   current chromosome: 18
#  #   current chromosome: 19
#  #   current chromosome: 20
#  #   current chromosome: 21
#  #   current chromosome: 22
#  #   current chromosome: X
#  #   current chromosome: Y
#  # Segments are written to: tumor01_cbs.seg
#  # Plotting

## ----eval=FALSE---------------------------------------------------------------
#  plotMosdepth_t(bed = "tumor.regions.bed.gz")

## -----------------------------------------------------------------------------
sessionInfo()

