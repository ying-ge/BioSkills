## ----global_options, include = FALSE------------------------------------------
library(knitr)
knitr::opts_chunk$set(fig.path='figs/', warning=FALSE, message=FALSE, collapse=TRUE)
source("render_toc.R")

## ----toc, echo = FALSE--------------------------------------------------------
render_toc("circRNAprofiler.Rmd")

## ----echo=FALSE, out.width='70%', fig.align="center", fig.cap="\\label{fig:figs} Figure 1: Schematic representation of the circRNA analysis workflow implemented by circRNAprofiler. The grey boxes represent the 15 modules with the main R-functions reported in italics. The different type of sequences that can be selected are depicted in the dashed box. BSJ, Back-Spliced Junction."----

knitr::include_graphics("./images/image1.png")

## ----eval = FALSE-------------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE)){
#    install.packages("BiocManager")
#  }
#  

## ----eval = FALSE-------------------------------------------------------------
#  BiocManager::install("circRNAprofiler")

## ----eval = FALSE-------------------------------------------------------------
#  # The following initializes usage of Bioc devel
#  BiocManager::install(version='devel')
#  
#  BiocManager::install("circRNAprofiler")

## -----------------------------------------------------------------------------
library(circRNAprofiler)

# Packages needed for the vignettes
library(ggpubr)
library(ggplot2)
library(VennDiagram)
library(gridExtra)

## ----echo=FALSE, out.width='55%', fig.align="center", fig.cap="\\label{fig:figs} Figure 2: Example of a project folder structure"----
knitr::include_graphics("./images/image2.png")

## ----eval = FALSE-------------------------------------------------------------
#  initCircRNAprofiler(projectFolderName = "projectCirc", detectionTools =
#                        "mapsplice")

## ----eval = FALSE-------------------------------------------------------------
#  initCircRNAprofiler(
#      projectFolderName = "projectCirc",
#      detectionTools = c("mapsplice", "nclscan", "circmarker")
#  )

## ----echo=FALSE---------------------------------------------------------------
experiment <-
    read.table(
        "experiment.txt",
        header = TRUE,
        stringsAsFactors = FALSE,
        sep = "\t"
    )
head(experiment)

## ----echo=FALSE---------------------------------------------------------------
motifs <-
    read.table(
        "motifs.txt",
        stringsAsFactors = FALSE,
        header = TRUE,
        sep = "\t"
    )
head(motifs)


## ----echo=FALSE---------------------------------------------------------------
traits <-
    read.table(
        "traits.txt",
        stringsAsFactors = FALSE,
        header = TRUE,
        sep = "\t"
    )
head(traits)

## ----echo=FALSE---------------------------------------------------------------
miRs <-
    read.table(
        "miRs.txt",
        header = TRUE,
        stringsAsFactors = FALSE,
        sep = "\t"
    )
head(miRs)

## ----echo=FALSE---------------------------------------------------------------
transcripts <-
    read.table(
        "transcripts.txt",
        header = TRUE,
        stringsAsFactors = FALSE,
        sep = "\t"
    )
head(transcripts)


## ----echo=FALSE---------------------------------------------------------------
circRNApredictions <-
    read.table(
        "circRNAs_test.txt",
        header = TRUE,
        stringsAsFactors = TRUE,
        sep = "\t"
    )
head(circRNApredictions)


## ----eval=FALSE---------------------------------------------------------------
#  # Path to experiment.txt
#  pathToExperiment <- system.file("extdata", "experiment.txt",
#                                  package ="circRNAprofiler")
#  
#  

## ----eval = FALSE-------------------------------------------------------------
#  # Set project folder projectCirc as your working directory and run:
#  check <- checkProjectFolder()
#  check

## -----------------------------------------------------------------------------
# Set project folder projectCirc as your working directory.
# Download gencode.V19.annotation.gtf from https://www.gencodegenes.org/ and 
# put it in the working directory, then run:
# gtf <- formatGTF("gencode.V19.annotation.gtf")

# For example purpose load a short version of the formatted gtf file generated
# using the command above.
data("gtf")
head(gtf)


## -----------------------------------------------------------------------------
# Set working directory to projectCirc which contains a short version of the raw
# files containing the detected circRNAs. The run: 
# backSplicedJunctions <- getBackSplicedJunctions(gtf)

# Alternatively, you can load the object containing the whole set of circRNAs detected
# in the heart generated running the command above.
data("backSplicedJunctions")
head(backSplicedJunctions)

## ----fig.align="center", fig.width = 10, fig.height = 3-----------------------
# Plot
p <- ggplot(backSplicedJunctions, aes(x = tool)) +
    geom_bar() +
    labs(title = "", x = "Detection tool", y = "No. of circRNAs") +
    theme_classic()

# Run getDetectionTools() to get the code corresponding to the circRNA
# detection tools.
dt <- getDetectionTools() %>%
    dplyr::filter( name %in% c("mapsplice","nclscan", "circmarker"))%>%
    gridExtra::tableGrob(rows=NULL)

# Merge plots
gridExtra::grid.arrange(p, dt, nrow=1)

## -----------------------------------------------------------------------------
# If you set projectCirc as your working directory, then run:
# mergedBSJunctions <- mergeBSJunctions(backSplicedJunctions, gtf)

# Alternatively, you can load the object containing the whole set of circRNAs 
# detected in the heart merged using the code above.
data("mergedBSJunctions")
head(mergedBSJunctions)


## ----fig.align = "center", fig.width = 10, fig.height = 4---------------------
# Plot
p <- ggplot(mergedBSJunctions, aes(x = tool)) +
    geom_bar() +
    labs(title = "", x = "Detection tool", y = "No. of circRNAs") +
    theme_classic()

gridExtra::grid.arrange(p, dt, nrow=1)


## -----------------------------------------------------------------------------
# If you set projectCirc as your working directory, then run:
filteredCirc <-
filterCirc(mergedBSJunctions, allSamples = FALSE, min = 5)


## ----fig.align="center", fig.width = 10, fig.height = 4-----------------------
# Plot
p <- ggplot(filteredCirc, aes(x = tool)) +
    geom_bar() +
    labs(title = "", x = "Detection tool", y = "No. of circRNAs") +
    theme_classic()

gridExtra::grid.arrange(p, dt, nrow=1)

## ----fig.align="center", fig.width = 5, fig.height = 4------------------------
# Plot using Venn diagram
cm <- filteredCirc[base::grep("cm", filteredCirc$tool), ]
ms <- filteredCirc[base::grep("ms", filteredCirc$tool), ]
ns <- filteredCirc[base::grep("ns", filteredCirc$tool), ]

p <- VennDiagram::draw.triple.venn(
    area1 = length(cm$id),
    area2 = length(ms$id),
    area3 = length(ns$id),
    n12 = length(intersect(cm$id, ms$id)),
    n23 = length(intersect(ms$id, ns$id)),
    n13 = length(intersect(cm$id, ns$id)),
    n123 = length(Reduce(
        intersect, list(cm$id, ms$id, ns$id)
    )),
    category = c("cm", "ms", "ns"),
    lty = "blank",
    fill = c("skyblue", "pink1", "mediumorchid")
)


## -----------------------------------------------------------------------------
# Compare condition B Vs A
# If you set projectCirc as your working directory, then run:
deseqResBvsA <-
    getDeseqRes(
        filteredCirc,
        condition = "A-B",
        fitType = "local",
        pAdjustMethod = "BH"
    )
head(deseqResBvsA)





## -----------------------------------------------------------------------------
# Compare condition C Vs A
deseqResCvsA <-
    getDeseqRes(
        filteredCirc,
        condition = "A-C",
        fitType = "local",
        pAdjustMethod = "BH"
    )
head(deseqResCvsA)

## ----fig.align="center", fig.height= 8, fig.width = 8-------------------------
# We set the xlim and ylim to the same values for both plots to make them
# comparable. Before setting the axis limits, you should visualize the 
# plots with the default values to be able to define the correct limits.
# An error might occur due to the log10 transformation of the padj values 
# (e.g. log10(0) = inf). In that case set setyLim = TRUE and specify the the 
# y limit manually.
p1 <-
    volcanoPlot(
        deseqResBvsA,
        log2FC = 1,
        padj = 0.05,
        title = "DCMs Vs. Con",
        setxLim = TRUE,
        xlim = c(-8 , 7.5),
        setyLim = TRUE,
        ylim = c(0 , 4),
        gene = FALSE
    )
p2 <-
    volcanoPlot(
        deseqResCvsA,
        log2FC = 1,
        padj = 0.05,
        title = "HCMs Vs. Con",
        setxLim = TRUE,
        xlim = c(-8 , 7.5),
        setyLim = TRUE,
        ylim = c(0 , 4),
        gene = FALSE
    )
ggarrange(p1, 
          p2, 
          ncol = 1, 
          nrow = 2)

## ----eval = FALSE-------------------------------------------------------------
#  # Compare condition B Vs A
#  edgerResBvsA <-
#      getEdgerRes(
#          filteredCirc,
#          condition = "A-B",
#          normMethod = "TMM",
#          pAdjustMethod = "BH"    )
#  head(edgerResBvsA)

## ----eval = FALSE-------------------------------------------------------------
#  # Compare condition C Vs A
#  edgerResCvsA <-
#      getEdgerRes(
#          filteredCirc,
#          condition = "A-C",
#          normMethod = "TMM",
#          pAdjustMethod = "BH"
#      )
#  head(edgerResCvsA)

## ----eval = FALSE-------------------------------------------------------------
#  liftedBSJcoords <- liftBSJcoords(filteredCirc, map = "hg19ToMm9",
#                                   annotationHubID = "AH14155")

## -----------------------------------------------------------------------------
# If you want to analysis specific transcripts report them in transcripts.txt (optional).
# If transcripts.txt is not present in your working directory specify pathToTranscripts.
# As default transcripts.txt is searched in the wd.

# As an example of the 1458 filtered circRNAs we annotate only the firt 30 
# circRNAs
annotatedBSJs <- annotateBSJs(filteredCirc[1:30,], gtf, isRandom = FALSE) 
head(annotatedBSJs)

## -----------------------------------------------------------------------------
# First find frequency of single exon circRNAs
f <-
    sum((annotatedBSJs$exNumUpBSE == 1 |
            annotatedBSJs$exNumDownBSE == 1) ,
        na.rm = TRUE) / (nrow(annotatedBSJs) * 2)

# Retrieve random back-spliced junctions
randomBSJunctions <-
    getRandomBSJunctions(gtf, n = nrow(annotatedBSJs), f = f, setSeed = 123)
head(randomBSJunctions)

## ----eval = FALSE-------------------------------------------------------------
#  annotatedRBSJs <- annotateBSJs(randomBSJunctions, gtf, isRandom = TRUE)

## ----fig.align="center", fig.width = 13, fig.height = 8, eval = FALSE---------
#  # annotatedBSJs act as foreground data set
#  # annotatedRBSJs act as background data set
#  
#  # Length of flanking introns
#  p1 <- plotLenIntrons(
#      annotatedBSJs,
#      annotatedRBSJs,
#      title = "Length flanking introns",
#      df1Name = "predicted",
#      df2Name = "random",
#      setyLim = TRUE,
#      ylim = c(0,7)
#  )
#  
#  # Length of back-splided exons
#  p2 <- plotLenBSEs(
#      annotatedBSJs,
#      annotatedRBSJs,
#      title = "Length back-splided exons",
#      df1Name = "predicted",
#      df2Name = "random",
#      setyLim = TRUE,
#      ylim = c(0,7)
#  )
#  
#  # No. of circRNAs produced from the host genes
#  p3 <-
#      plotHostGenes(annotatedBSJs, title = "# CircRNAs produced from host genes")
#  
#  # No. of exons in between the back-spliced junctions
#  p4 <-
#      plotExBetweenBSEs(annotatedBSJs, title = "# Exons between back-spliced junctions")
#  
#  # Position of back-spliced exons within the host transcripts
#  p5 <-
#      plotExPosition(annotatedBSJs,
#          n = 1,
#          title = "Position back-spliced exons in the transcripts")
#  
#  # Total no. of exons within the host transcripts
#  p6 <-
#      plotTotExons(annotatedBSJs, title = " Total number of exons in the host transcripts")
#  
#  # Combine plots
#  ggarrange(p1,
#      p2,
#      p3,
#      p4,
#      p5,
#      p6,
#      ncol = 2,
#      nrow = 3)
#  

## ----echo=FALSE, out.width='100%', fig.align="center", fig.cap="\\label{fig:figs} Comparison of structural features extracted from the subset of 1458 filtered  back-spliced junctions compared to an equal number of randomly generated back-spliced junctions."----
knitr::include_graphics("./images/image3.png")

## ----eval = FALSE-------------------------------------------------------------
#  # Select ALPK2:-:chr18:56247780:56246046 circRNA
#  annotatedCirc <-
#  annotatedBSJs[annotatedBSJs$id == "ALPK2:-:chr18:56247780:56246046", ]
#  
#  # As background data set we used all the remaining 1457 filered circRNAs.
#  # Alternatively the subset of randomly generated back-spliced junctions can be used.
#  annotatedBackgroundCircs <-
#  annotatedBSJs[which(annotatedBSJs$id != "ALPK2:-:chr18:56247780:56246046"), ]

## -----------------------------------------------------------------------------
# All the sequences will be retrieved from the BSgenome package which contains 
# the serquences of the genome of interest
if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)){
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
}

# Get genome
genome <- BSgenome::getBSgenome("BSgenome.Hsapiens.UCSC.hg19")

## ----eval = FALSE-------------------------------------------------------------
#  
#  # Foreground target sequences
#  targetsFTS_circ <-
#  getCircSeqs(annotatedCirc, gtf, genome)
#  
#  # Background target sequences
#  targetsBTS_circ <-
#  getCircSeqs(annotatedBackgroundCircs, gtf, genome)
#  

## ----eval = FALSE-------------------------------------------------------------
#  # Foreground target sequences
#  targetsFTS_bsj <-
#      getSeqsAcrossBSJs(annotatedCirc, gtf, genome)
#  
#  # Background target sequences
#  targetsBTS_bsj <-
#  getCircSeqs(annotatedBackgroundCircs, gtf, genome)
#  

## ----eval = FALSE-------------------------------------------------------------
#  # Foreground target sequences
#  targetsFTS_gr <-
#      getSeqsFromGRs(
#          annotatedCirc,
#          genome,
#          lIntron = 200,
#          lExon = 9,
#          type = "ie"
#          )
#  # Background target sequences.
#  targetsBTS_gr <-
#      getSeqsFromGRs(
#          annotatedBackgroundCircs,
#          genome,
#          lIntron = 200,
#          lExon = 9,
#          type = "ie")

## ----eval = FALSE-------------------------------------------------------------
#  # If you want to analysis also custom motifs report them in motifs.txt (optional).
#  # If motifs.txt is not present in your working directory specify pathToMotifs.
#  # As default motifs.txt is searched in the wd.
#  
#  # E.g. in motifs.txt we added RBM20 consensus motif since it is not present in
#  # the ATtRACT database.
#  
#  # Find motifs in the foreground target sequences
#  motifsFTS_gr <-
#      getMotifs(targetsFTS_gr,
#                width = 6,
#                database = 'ATtRACT',
#                species = "Hsapiens",
#                rbp = TRUE,
#                reverse = FALSE)
#  # Find motifs in the background target sequences
#  motifsBTS_gr <-
#      getMotifs(targetsBTS_gr,
#                width = 6,
#                database = 'ATtRACT',
#                species = "Hsapiens",
#                rbp = TRUE,
#                reverse = FALSE)

## ----eval = FALSE-------------------------------------------------------------
#  mergedMotifsFTS_gr <- mergeMotifs(motifsFTS_gr)
#  mergedMotifsBTS_gr <- mergeMotifs(motifsBTS_gr)

## ----fig.align="center", fig.width = 7, fig.height = 7, eval = FALSE----------
#  # Plot log2FC and normalized counts
#  
#  # Normalize using the number of target sequences. You can do this normalization if you
#  # analyzed the same number of target sequences and extracted the same number of
#  # nucleotides from the latter.
#  
#  # Alternatively you can use the length of the target sequences.
#  # If you are analyzing the flanking upstream and downstream sequences then normalize
#  # using the length of the latters.E.g.:
#  # nf1 = sum(targetsFTS_gr$upGR$length, na.rm = T)+sum(targetsFTS_gr$downGR$length, na.rm = T)
#  # nf2 = sum(targetsBTS_gr$upGR$length, na.rm = T)+sum(targetsBTS_gr$downGR$length, na.rm = T)
#  
#  # If you are analyzing the circRNA sequences then normalize using the length of
#  # the latter. E.g.:
#  # nf1 = sum(targetsFTS_circ$circ$length, na.rm = T)
#  # nf2 = sum(targetsBTS_circ$circ$length, na.rm = T)
#  
#  # If you are analyzing the BSJ then normalize using the length of the latter. E.g.:
#  # nf1 = sum(targetsFTS_bsj$bsj$length, na.rm = T)
#  # nf2 = sum(targetsBTS_bsj$bsj$length, na.rm = T)
#  
#  p <-
#      plotMotifs(
#          mergedMotifsFTS_gr,
#          mergedMotifsBTS_gr,
#          nf1 = sum(targetsFTS_gr$upGR$length, na.rm = T)+sum(targetsFTS_gr$upGR$length, na.rm = T) ,
#          nf2 = sum(targetsBTS_gr$upGR$length, na.rm = T)+sum(targetsBTS_gr$upGR$length, na.rm = T),
#          log2FC = 1,
#          removeNegLog2FC = TRUE,
#          df1Name = "circALPK2",
#          df2Name = "Other circRNAs",
#          angle = 45
#      )
#  ggarrange(p[[1]],
#            p[[2]],
#            labels = c("", ""),
#            ncol = 2,
#            nrow = 1)
#  

## ----echo=FALSE, out.width='70%', fig.align="center", fig.cap="\\label{fig:figs} Bar chart showing the log2FC (cut-off = 1) and the normalized counts of the RBP motifs found in the region flanking the predicted back-spliced junction of circALPK2 compared to the remaining subset of 1457 filtered circRNAs."----
knitr::include_graphics("./images/image4.png")

## ----eval = FALSE-------------------------------------------------------------
#  # Type p[[3]] to get the table
#  head(p[[3]])

## ----eval = FALSE-------------------------------------------------------------
#  # If you want to analysis only a subset of miRs, then specify the miR ids in
#  # miRs.txt (optional).
#  # If miRs.txt is not present in your working directory specify pathToMiRs.
#  # As default miRs.txt is searched in the wd.
#  miRsites <-
#      getMiRsites(
#          targetsFTS_circ,
#          miRspeciesCode = "hsa",
#          miRBaseLatestRelease = TRUE,
#          totalMatches = 6,
#          maxNonCanonicalMatches = 1
#      )

## ----eval = FALSE-------------------------------------------------------------
#  rearragedMiRres <- rearrangeMiRres(miRsites)

## ----eval=FALSE---------------------------------------------------------------
#  # If multiple circRNAs have been analyzed for the presence of miR binding sites
#  # the following code can store the predictions for each circRNA in a
#  # different sheet of an xlsx file for a better user consultation.
#  i <- 1
#  j <- 1
#  while (i <= (length(rearragedMiRres))) {
#      write.xlsx2(
#          rearragedMiRres[[i]][[1]],
#          "miRsites_TM6_NCM1.xlsx",
#          paste("sheet", j, sep = ""),
#          append = TRUE
#      )
#      j <- j + 1
#      write.xlsx2(
#          rearragedMiRres[[i]][[2]],
#          "miRsites_TM6_NCM1.xlsx",
#          paste("sheet", j, sep = ""),
#          append = TRUE
#      )
#      i <- i + 1
#      j <- j + 1
#  }
#  

## ----eval = FALSE-------------------------------------------------------------
#  # Plot miRNA analysis results
#  
#  p <- plotMiR(rearragedMiRres,
#               n = 40,
#               color = "blue",
#               miRid = TRUE,
#               id = 1)
#  p

## ----eval = FALSE-------------------------------------------------------------
#  snpsGWAS <-
#      annotateSNPsGWAS(targetsFTS_gr, assembly = "hg19", makeCurrent = TRUE)
#  

## ----eval = FALSE-------------------------------------------------------------
#  repeats <-
#      annotateRepeats(targetsFTS_gr, annotationHubID = "AH5122", complementary = TRUE)
#  

## -----------------------------------------------------------------------------
sessionInfo()

