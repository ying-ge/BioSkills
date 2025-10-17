## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
set.seed(2)

## -----------------------------------------------------------------------------
library(glmGamPoi)

## -----------------------------------------------------------------------------
# overdispersion = 1/size
counts <- rnbinom(n = 10, mu = 5, size = 1/0.7)

# design = ~ 1 means that an intercept-only model is fit
fit <- glm_gp(counts, design = ~ 1)
fit

# Internally fit is just a list:
as.list(fit)[1:2]

## ----warning=FALSE, message = FALSE-------------------------------------------
library(SummarizedExperiment)
library(DelayedMatrixStats)

## -----------------------------------------------------------------------------
# The full dataset with 33,000 genes and 4340 cells
# The first time this is run, it will download the data
pbmcs <- TENxPBMCData::TENxPBMCData("pbmc4k")

# I want genes where at least some counts are non-zero
non_empty_rows <- which(rowSums2(assay(pbmcs)) > 0)
pbmcs_subset <- pbmcs[sample(non_empty_rows, 300), ]
pbmcs_subset

## -----------------------------------------------------------------------------
fit <- glm_gp(pbmcs_subset, on_disk = FALSE)
summary(fit)

## ----warning=FALSE------------------------------------------------------------
# Explicitly realize count matrix in memory so that it is a fair comparison
pbmcs_subset <- as.matrix(assay(pbmcs_subset))
model_matrix <- matrix(1, nrow = ncol(pbmcs_subset))


bench::mark(
  glmGamPoi_in_memory = {
    glm_gp(pbmcs_subset, design = model_matrix, on_disk = FALSE)
  }, glmGamPoi_on_disk = {
    glm_gp(pbmcs_subset, design = model_matrix, on_disk = TRUE)
  }, DESeq2 = suppressMessages({
    dds <- DESeq2::DESeqDataSetFromMatrix(pbmcs_subset,
                        colData = data.frame(name = seq_len(4340)),
                        design = ~ 1)
    dds <- DESeq2::estimateSizeFactors(dds, "poscounts")
    dds <- DESeq2::estimateDispersions(dds, quiet = TRUE)
    dds <- DESeq2::nbinomWaldTest(dds, minmu = 1e-6)
  }), edgeR = {
    edgeR_data <- edgeR::DGEList(pbmcs_subset)
    edgeR_data <- edgeR::calcNormFactors(edgeR_data)
    edgeR_data <- edgeR::estimateDisp(edgeR_data, model_matrix)
    edgeR_fit <- edgeR::glmFit(edgeR_data, design = model_matrix)
  }, check = FALSE, min_iterations = 3
)

## ----message=FALSE, warning=FALSE---------------------------------------------
# Results with my method
fit <- glm_gp(pbmcs_subset, design = model_matrix, on_disk = FALSE)

# DESeq2
dds <- DESeq2::DESeqDataSetFromMatrix(pbmcs_subset, 
                        colData = data.frame(name = seq_len(4340)),
                        design = ~ 1)
sizeFactors(dds)  <- fit$size_factors
dds <- DESeq2::estimateDispersions(dds, quiet = TRUE)
dds <- DESeq2::nbinomWaldTest(dds, minmu = 1e-6)

#edgeR
edgeR_data <- edgeR::DGEList(pbmcs_subset, lib.size = fit$size_factors)
edgeR_data <- edgeR::estimateDisp(edgeR_data, model_matrix)
edgeR_fit <- edgeR::glmFit(edgeR_data, design = model_matrix)

## ----coefficientComparison, fig.height=5, fig.width=10, warning=FALSE, echo = FALSE----
par(mfrow = c(2, 4), cex.main = 2, cex.lab = 1.5)
plot(fit$Beta[,1], coef(dds)[,1] / log2(exp(1)), pch = 16, 
     main = "Beta Coefficients", xlab = "glmGamPoi", ylab = "DESeq2")
abline(0,1)
plot(fit$Beta[,1], edgeR_fit$unshrunk.coefficients[,1], pch = 16,
     main = "Beta Coefficients", xlab = "glmGamPoi", ylab = "edgeR")
abline(0,1)

plot(fit$Mu[,1], assay(dds, "mu")[,1], pch = 16, log="xy",
     main = "Gene Mean", xlab = "glmGamPoi", ylab = "DESeq2")
abline(0,1)
plot(fit$Mu[,1], edgeR_fit$fitted.values[,1], pch = 16, log="xy",
     main = "Gene Mean", xlab = "glmGamPoi", ylab = "edgeR")
abline(0,1)

plot(fit$overdispersions, rowData(dds)$dispGeneEst, pch = 16, log="xy",
     main = "Overdispersion", xlab = "glmGamPoi", ylab = "DESeq2")
abline(0,1)
plot(fit$overdispersions, edgeR_fit$dispersion, pch = 16, log="xy",
     main = "Overdispersion", xlab = "glmGamPoi", ylab = "edgeR")
abline(0,1)


## -----------------------------------------------------------------------------
sce <- muscData::Kang18_8vs8()
colData(sce)

## -----------------------------------------------------------------------------
set.seed(1)
# Take highly expressed genes and proper cells:
sce_subset <- sce[rowSums(counts(sce)) > 100, 
                  sample(which(sce$multiplets == "singlet" & 
                              ! is.na(sce$cell) &
                              sce$cell %in% c("CD4 T cells", "B cells", "NK cells")), 
                         1000)]
# Convert counts to dense matrix
counts(sce_subset) <- as.matrix(counts(sce_subset))
# Remove empty levels because glm_gp() will complain otherwise
sce_subset$cell <- droplevels(sce_subset$cell)

## -----------------------------------------------------------------------------
fit <- glm_gp(sce_subset, design = ~ cell + stim +  stim:cell - 1,
              reference_level = "NK cells")
summary(fit)

## -----------------------------------------------------------------------------
colnames(fit$Beta)

## -----------------------------------------------------------------------------
# The contrast argument specifies what we want to compare
# We test the expression difference of stimulated and control T-cells
#
# There is no sample label in the colData, so we create it on the fly
# from `stim` and `ind` columns in colData(fit$data).
de_res <- test_de(fit, contrast = `stimstim` + `cellCD4 T cells:stimstim`, 
                  pseudobulk_by = paste0(stim, "-", ind)) 

# The large `lfc` values come from groups were nearly all counts are 0
# Setting them to Inf makes the plots look nicer
de_res$lfc <- ifelse(abs(de_res$lfc) > 20, sign(de_res$lfc) * Inf, de_res$lfc)

# Most different genes
head(de_res[order(de_res$pval), ])

## -----------------------------------------------------------------------------
library(ggplot2)
ggplot(de_res, aes(x = lfc, y = -log10(pval))) +
  geom_point(size = 0.6, aes(color = adj_pval < 0.1)) +
  ggtitle("Volcano Plot", "Genes that change most through interferon-beta treatment in T cells")

## -----------------------------------------------------------------------------
marker_genes <- test_de(fit, `cellCD4 T cells` - `cellB cells`, sort_by = pval)
head(marker_genes)

## -----------------------------------------------------------------------------
marker_genes2 <- test_de(fit, (`cellCD4 T cells` + `cellCD4 T cells:stimstim`) - 
                               (`cellB cells` + `cellB cells:stimstim`), 
                        sort_by = pval)

head(marker_genes2)

## -----------------------------------------------------------------------------
# Create a data.frame with the expression values, gene names, and cell types
tmp <- data.frame(gene = rep(marker_genes$name[1:6], times = ncol(sce_subset)),
                  expression = c(counts(sce_subset)[marker_genes$name[1:6], ]),
                  celltype = rep(sce_subset$cell, each = 6))

ggplot(tmp, aes(x = celltype, y = expression)) +
  geom_jitter(height = 0.1) +
  stat_summary(geom = "crossbar", fun = "mean", color = "red") +
  facet_wrap(~ gene, scales = "free_y") +
  ggtitle("Marker genes of B vs. T cells")

## -----------------------------------------------------------------------------
sessionInfo()

