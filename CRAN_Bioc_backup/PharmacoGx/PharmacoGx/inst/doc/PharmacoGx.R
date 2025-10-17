## ----setup, include=FALSE, cache=FALSE, message = FALSE-----------------------

library("knitr")

### Chunk options: see http://yihui.name/knitr/options/ ###

## Text results
opts_chunk$set(echo = TRUE, warning = TRUE, message = TRUE, include = TRUE)

## Code decoration
opts_chunk$set(tidy = FALSE, comment = NA, highlight = TRUE)


## ----knitcitation, include=FALSE----------------------------------------------
library(knitcitations)
cleanbib()
cite_options(citation_format = "pandoc")

## ----install-pkg, eval=FALSE, results='hide'----------------------------------
#  if (!requireNamespace("BiocManager", quietly=TRUE))
#      install.packages("BiocManager")
#  BiocManager::install('PharmacoGx')

## ----eval=FALSE---------------------------------------------------------------
#  library(PharmacoGx)

## ----download-example, eval=FALSE---------------------------------------------
#  availablePSets()
#  GDSC <- downloadPSet("GDSC_2020(v2-8.2)")

## ----inconsistencies, results='hide', eval=TRUE, message = FALSE--------------
  library(Biobase)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(PharmacoGx)
  data("GDSCsmall")
  data("CCLEsmall")
  commonGenes <- intersect(fNames(GDSCsmall, "rna"),
                           fNames(CCLEsmall,"rna"))
  common <- intersectPSet(list('CCLE'=CCLEsmall,
                               'GDSC'=GDSCsmall),
                          intersectOn=c("cell.lines", "drugs"),
                          strictIntersect=TRUE)


  GDSC.auc <- summarizeSensitivityProfiles(
                common$GDSC,
                sensitivity.measure='auc_published',
                summary.stat="median",
                verbose=FALSE)
  CCLE.auc <- summarizeSensitivityProfiles(
                common$CCLE,
                sensitivity.measure='auc_published',
                summary.stat="median",
                verbose=FALSE)

  GDSC.ic50 <- summarizeSensitivityProfiles(
                common$GDSC,
                sensitivity.measure='ic50_published',
                summary.stat="median",
                verbose=FALSE)
  CCLE.ic50 <- summarizeSensitivityProfiles(
                common$CCLE,
                sensitivity.measure='ic50_published',
                summary.stat="median",
                verbose=FALSE)

  GDSCexpression <- summarizeMolecularProfiles(common$GDSC,
                                        cellNames(common$GDSC),
                                        mDataType="rna",
                                        features=commonGenes,
                                        verbose=FALSE)
  CCLEexpression <- summarizeMolecularProfiles(common$CCLE,
                                         cellNames(common$CCLE),
                                         mDataType="rna",
                                         features=commonGenes,
                                         verbose=FALSE)
  gg <- fNames(common[[1]], 'rna')
  cc <- cellNames(common[[1]])

  ge.cor <- sapply(cc, function (x, d1, d2) {
    stats::cor(d1[ , x], d2[ , x], method="spearman",
                use="pairwise.complete.obs")
  ## TO DO:: Ensure all assays are name so we can call by name instead of index
  }, d1=assay(GDSCexpression, 1), d2=assay(CCLEexpression, 1))
  ic50.cor <- sapply(cc, function (x, d1, d2) {
    stats::cor(d1[, x], d2[ , x], method="spearman",
                use="pairwise.complete.obs")
  }, d1=GDSC.ic50, d2=CCLE.ic50)
  auc.cor <- sapply(cc, function (x, d1, d2) {
    stats::cor(d1[ , x], d2[ , x], method="spearman",
                use="pairwise.complete.obs")
  }, d1=GDSC.auc, d2=CCLE.auc)

  w1 <- stats::wilcox.test(x=ge.cor, y=auc.cor,
                           conf.int=TRUE, exact=FALSE)
  w2 <- stats::wilcox.test(x=ge.cor, y=ic50.cor,
                           conf.int=TRUE, exact=FALSE)
  yylim <- c(-1, 1)
  ss <- sprintf("GE vs. AUC = %.1E\nGE vs. IC50 = %.1E",
                w1$p.value, w2$p.value)

## ----fig1---------------------------------------------------------------------
  boxplot(list("GE"=ge.cor,
               "AUC"=auc.cor,
               "IC50"=ic50.cor),
          main="Concordance between cell lines",
          ylab=expression(R[s]),
          sub=ss,
          ylim=yylim,
          col="lightgrey",
          pch=20,
          border="black")

## ----eval=TRUE, results='asis'------------------------------------------------
  library(PharmacoGx)
  library(pander)
  data(CMAPsmall)
  drug.perturbation <- drugPerturbationSig(CMAPsmall,
                                           mDataType="rna",
                                           verbose=FALSE)
  data(HDAC_genes)

  res <- apply(drug.perturbation[,,c("tstat", "fdr")],
               2, function(x, HDAC){
	    return(connectivityScore(x=x,
	                             y=HDAC[,2,drop=FALSE],
                               method="fgsea", nperm=100))
	}, HDAC=HDAC_genes)

  rownames(res) <- c("Connectivity", "P Value")
  res <- t(res)
  res <- res[order(res[,1], decreasing=TRUE),]
  pander::pandoc.table(res,
    caption='Connectivity Score results for HDAC inhibitor gene signature.',
    style = 'rmarkdown')

## ----biomarkers, eval=TRUE, results='asis'------------------------------------
  library(pander)
  data(CCLEsmall)
  features <- fNames(CCLEsmall, "rna")[
                          which(featureInfo(CCLEsmall,
                                            "rna")$Symbol == "NQO1")]
  sig.rna <- drugSensitivitySig(object=CCLEsmall,
                            mDataType="rna",
                            drugs=c("17-AAG"),
                            features=features,
                            sensitivity.measure="auc_published",
                            molecular.summary.stat="median",
                            sensitivity.summary.stat="median",
                            verbose=FALSE)
  sig.mut <- drugSensitivitySig(object=CCLEsmall,
                            mDataType="mutation",
                            drugs=c("PD-0325901"),
                            features="BRAF",
                            sensitivity.measure="auc_published",
                            molecular.summary.stat="and",
                            sensitivity.summary.stat="median",
                            verbose=FALSE)
  sig <- rbind(sig.rna, sig.mut)
  rownames(sig) <- c("17-AAG + NQO1","PD-0325901 + BRAF")
  colnames(sig) <- dimnames(sig.mut)[[3]]
  pander::pandoc.table(t(sig), style = "rmarkdown", caption='P Value of Gene-Drug Association' )

## ----sessioninfo, eval = TRUE-------------------------------------------------
# set eval = FALSE if you don't want this info (useful for reproducibility) to appear
sessionInfo()

