## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>",
    message = FALSE,
    warning = FALSE,
    fig.align = "center",
    crop = NULL
)
library(BiocStyle)

## ----install-bioc,eval=FALSE--------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE)) {
#      install.packages("BiocManager")
#  }
#  
#  BiocManager::install("microbiomeMarker")

## ----install-gh,eval=FALSE----------------------------------------------------
#  if (!requireNamespace("remotes", quietly = TRUE)) {
#      install.packages("remotes")
#  }
#  remotes::install_github("yiluheihei/microbiomeMarker")

## ----load---------------------------------------------------------------------
library(microbiomeMarker)

## ----import-dada2-------------------------------------------------------------
seq_tab <- readRDS(
    system.file(
        "extdata", "dada2_seqtab.rds",
        package = "microbiomeMarker"
    )
)
tax_tab <- readRDS(
    system.file(
        "extdata", "dada2_taxtab.rds",
        package = "microbiomeMarker"
    )
)
sam_tab <- read.table(
    system.file(
        "extdata", "dada2_samdata.txt",
        package = "microbiomeMarker"
    ),
    sep = "\t",
    header = TRUE,
    row.names = 1
)
ps <- import_dada2(seq_tab = seq_tab, tax_tab = tax_tab, sam_tab = sam_tab)
ps

## ----import-qiime2,message=FALSE----------------------------------------------
otuqza_file <- system.file(
    "extdata", "table.qza",
    package = "microbiomeMarker"
)
taxaqza_file <- system.file(
    "extdata", "taxonomy.qza",
    package = "microbiomeMarker"
)
sample_file <- system.file(
    "extdata", "sample-metadata.tsv",
    package = "microbiomeMarker"
)
treeqza_file <- system.file(
    "extdata", "tree.qza",
    package = "microbiomeMarker"
)

ps <- import_qiime2(
    otu_qza = otuqza_file, taxa_qza = taxaqza_file,
    sam_tab = sample_file, tree_qza = treeqza_file
)
ps

## ----norm---------------------------------------------------------------------
# take tss as example
norm_tss(ps)

normalize(ps, method = "TSS")

## ----norm-note,eval=FALSE-----------------------------------------------------
#  data(kostic_crc)
#  mm_test <- normalize(kostic_crc, method = "CPM") %>%
#      run_lefse(
#          wilcoxon_cutoff = 0.01,
#          norm = "none", # must be "none" since the input has been normalized
#          group = "DIAGNOSIS",
#          kw_cutoff = 0.01,
#          multigrp_strat = TRUE,
#          lda_cutoff = 4
#      )
#  # equivalent to
#  run_lefse(
#      wilcoxon_cutoff = 0.01,
#      norm = "CPM",
#      group = "DIAGNOSIS",
#      kw_cutoff = 0.01,
#      multigrp_strat = TRUE,
#      lda_cutoff = 4
#  )

## ----two-group-test-----------------------------------------------------------
data(enterotypes_arumugam)
tg_welch <- run_test_two_groups(
    enterotypes_arumugam,
    group = "Gender",
    method = "welch.test"
)

# three significantly differential genera (marker)
tg_welch

# details of result of the three markers
head(marker_table(tg_welch))

## ----multi-group-test---------------------------------------------------------
# three groups
ps <- phyloseq::subset_samples(
    enterotypes_arumugam,
    Enterotype %in% c("Enterotype 3", "Enterotype 2", "Enterotype 1")
)
mg_anova <- run_test_multiple_groups(
    ps,
    group = "Enterotype",
    method = "anova"
)

# 24 markers
mg_anova

head(marker_table(mg_anova))

## ----edger--------------------------------------------------------------------
# contrast must be specified for two groups comparison
data(pediatric_ibd)
mm_edger <- run_edger(
    pediatric_ibd,
    group = "Class",
    pvalue_cutoff = 0.1,
    p_adjust = "fdr"
)
mm_edger

# multiple groups
data(cid_ying)
cid <- phyloseq::subset_samples(
    cid_ying,
    Consistency %in% c("formed stool", "liquid", "semi-formed")
)
mm_edger_mg <- run_edger(
    cid,
    group = "Consistency",
    method  = "QLFT",
    pvalue_cutoff = 0.05,
    p_adjust = "fdr"
)
mm_edger_mg

## ----lefse--------------------------------------------------------------------
data(kostic_crc)
kostic_crc_small <- phyloseq::subset_taxa(
    kostic_crc,
    Phylum %in% c("Firmicutes")
)
mm_lefse <- run_lefse(
    kostic_crc_small,
    wilcoxon_cutoff = 0.01,
    group = "DIAGNOSIS",
    kw_cutoff = 0.01,
    multigrp_strat = TRUE,
    lda_cutoff = 4
)

mm_lefse
head(marker_table(mm_lefse))

## ----rf-----------------------------------------------------------------------
# must specify the importance para for random forest
set.seed(2021)
# small example phyloseq object for test
ps_small <- phyloseq::subset_taxa(
    enterotypes_arumugam,
    Phylum %in% c("Firmicutes", "Bacteroidetes")
)
mm_lr <- run_sl(
    ps_small,
    group = "Gender",
    nfolds = 2,
    nrepeats = 1,
    taxa_rank = "Genus",
    top_n = 15,
    norm = "TSS",
    method = "LR",
)

marker_table(mm_lr)

## ----post-hoc-test------------------------------------------------------------
pht <- run_posthoc_test(ps, group = "Enterotype")
pht

# 24 significantly differential genera
markers <- marker_table(mg_anova)$feature
markers

# take a marker "p__Bacteroidetes|g__Bacteroides"
# for example, we will show "p__Bacteroidetes|g__Bacteroides"  differ from
# between Enterotype 2-Enterotype 1 and Enterotype 3-Enterotype 2.
extract_posthoc_res(pht, "p__Bacteroidetes|g__Bacteroides")[[1]]

## ----pair-wise-linear---------------------------------------------------------
# comparison between Enterotype 3 and Enterotype 2
mm_lv_pair <- run_limma_voom(
    ps,
    "Enterotype",
    contrast = c("Enterotype 3", "Enterotype 2"),
    pvalue_cutoff = 0.05,
    p_adjust = "fdr"
)
mm_lv_pair
head(marker_table(mm_lv_pair))

## ----plot-abundance-----------------------------------------------------------
p_abd <- plot_abundance(mm_lefse, group = "DIAGNOSIS")
p_abd

# customize the plot with ggplot2, modify the fill color manually
library(ggplot2)
p_abd + scale_fill_manual(values = c("Healthy" = "grey", "Tumor" = "red"))

## ----heatmap------------------------------------------------------------------
plot_heatmap(mm_edger, transform = "log10p", group = "Class")

## ----ef-plot------------------------------------------------------------------
# bar plot
plot_ef_bar(mm_lefse)

# dot plot
plot_ef_dot(mm_lefse)

## ----ef-plot-diff-------------------------------------------------------------
# set the x axis to log2 Fold Change automatically without manual intervention
plot_ef_bar(mm_edger)

## ----cladogram,fig.width=7,fig.height=7---------------------------------------
plot_cladogram(mm_lefse, color = c(Healthy = "darkgreen", Tumor = "red")) +
    theme(plot.margin = margin(0, 0, 0, 0))

## ----auc-roc------------------------------------------------------------------
set.seed(2021)
plot_sl_roc(mm_lr, group = "Gender")

## ----plot-pht-----------------------------------------------------------------
p_pht <- plot_postHocTest(pht, feature = "p__Bacteroidetes|g__Bacteroides")
p_pht

## ----customize-p-pht----------------------------------------------------------
p_pht & theme_bw()

## -----------------------------------------------------------------------------
sessionInfo()

