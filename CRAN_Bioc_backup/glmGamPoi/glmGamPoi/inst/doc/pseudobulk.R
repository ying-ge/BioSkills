## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
set.seed(2)

## ----out.width = '70%', echo = FALSE------------------------------------------
knitr::include_graphics("img/kang_data_overview.png")

## -----------------------------------------------------------------------------
library(SummarizedExperiment)
library(SingleCellExperiment)
sce <- muscData::Kang18_8vs8() 
# Keep only genes with more than 5 counts
sce <- sce[rowSums(counts(sce)) > 5,]
colData(sce)

## -----------------------------------------------------------------------------
library(glmGamPoi)
reduced_sce <- pseudobulk(sce, group_by = vars(ind, condition = stim, cell), 
                          fraction_singlet = mean(multiplets == "singlet"), n_cells = n())
colData(reduced_sce)

## -----------------------------------------------------------------------------
library(dplyr, warn.conflicts = FALSE)
colData(sce) %>%
  as_tibble() %>%
  group_by(ind, condition = stim, cell) %>%
  summarize(n_cells = n(), .groups = "drop") 

## -----------------------------------------------------------------------------
# Remove NA's
reduced_sce <- reduced_sce[,!is.na(reduced_sce$cell)]
# Use DESeq2's size factor calculation procedure
fit <- glm_gp(reduced_sce, design = ~ condition*cell + ind, size_factor = "ratio", verbose = TRUE)
res <- test_de(fit, contrast = cond(cell = "B cells", condition = "stim") - cond(cell = "B cells", condition = "ctrl"))

## -----------------------------------------------------------------------------
library(ggplot2, warn.conflicts = FALSE)
ggplot(res, aes(x = lfc, y = - log10(pval))) +
  geom_point(aes(color = adj_pval < 0.01), size = 0.5)

