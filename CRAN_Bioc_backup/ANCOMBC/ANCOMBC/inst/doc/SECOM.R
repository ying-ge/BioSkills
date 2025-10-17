## ----setup, message = FALSE, warning = FALSE, comment = NA--------------------
knitr::opts_chunk$set(message = FALSE, warning = FALSE, comment = NA, 
                      fig.width = 6.25, fig.height = 5)
library(ANCOMBC)
library(tidyverse)

## ----helper-------------------------------------------------------------------
get_upper_tri = function(cormat){
    cormat[lower.tri(cormat)] = NA
    diag(cormat) = NA
    return(cormat)
}

## ----getPackage, eval=FALSE---------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("ANCOMBC")

## ----load, eval=FALSE---------------------------------------------------------
#  library(ANCOMBC)

## -----------------------------------------------------------------------------
data(atlas1006, package = "microbiome")
tse = mia::makeTreeSummarizedExperimentFromPhyloseq(atlas1006)

# subset to baseline
tse = tse[, tse$time == 0]

# Re-code the bmi group
tse$bmi = recode(tse$bmi_group,
                 obese = "obese",
                 severeobese = "obese",
                 morbidobese = "obese")
# Subset to lean, overweight, and obese subjects
tse = tse[, tse$bmi %in% c("lean", "overweight", "obese")]

# Create the region variable
tse$region = recode(as.character(tse$nationality),
                    Scandinavia = "NE", UKIE = "NE", SouthEurope = "SE", 
                    CentralEurope = "CE", EasternEurope = "EE",
                    .missing = "unknown")

# Discard "EE" as it contains only 1 subject
# Discard subjects with missing values of region
tse = tse[, ! tse$region %in% c("EE", "unknown")]

print(tse)

## -----------------------------------------------------------------------------
set.seed(123)
# Linear relationships
res_linear = secom_linear(data = list(tse), assay_name = "counts",
                          tax_level = "Phylum", pseudo = 0, 
                          prv_cut = 0.5, lib_cut = 1000, corr_cut = 0.5, 
                          wins_quant = c(0.05, 0.95), method = "pearson", 
                          soft = FALSE, thresh_len = 20, n_cv = 10, 
                          thresh_hard = 0.3, max_p = 0.005, n_cl = 2)

# Nonlinear relationships
res_dist = secom_dist(data = list(tse), assay_name = "counts",
                      tax_level = "Phylum", pseudo = 0, 
                      prv_cut = 0.5, lib_cut = 1000, corr_cut = 0.5, 
                      wins_quant = c(0.05, 0.95), R = 1000, 
                      thresh_hard = 0.3, max_p = 0.005, n_cl = 2)

## -----------------------------------------------------------------------------
corr_linear = res_linear$corr_th
cooccur_linear = res_linear$mat_cooccur

# Filter by co-occurrence
overlap = 10
corr_linear[cooccur_linear < overlap] = 0

df_linear = data.frame(get_upper_tri(corr_linear)) %>%
  rownames_to_column("var1") %>%
  pivot_longer(cols = -var1, names_to = "var2", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(value = round(value, 2))

tax_name = sort(union(df_linear$var1, df_linear$var2))
df_linear$var1 = factor(df_linear$var1, levels = tax_name)
df_linear$var2 = factor(df_linear$var2, levels = tax_name)

heat_linear_th = df_linear %>%
  ggplot(aes(var2, var1, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", na.value = "grey",
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name = NULL) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  geom_text(aes(var2, var1, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Pearson (Thresholding)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1, 
                                   face = "italic"),
        axis.text.y = element_text(size = 12, face = "italic"),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 15),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") +
  coord_fixed()

heat_linear_th

## -----------------------------------------------------------------------------
corr_linear = res_linear$corr_fl
cooccur_linear = res_linear$mat_cooccur

# Filter by co-occurrence
overlap = 10
corr_linear[cooccur_linear < overlap] = 0

df_linear = data.frame(get_upper_tri(corr_linear)) %>%
  rownames_to_column("var1") %>%
  pivot_longer(cols = -var1, names_to = "var2", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(value = round(value, 2))

tax_name = sort(union(df_linear$var1, df_linear$var2))
df_linear$var1 = factor(df_linear$var1, levels = tax_name)
df_linear$var2 = factor(df_linear$var2, levels = tax_name)

heat_linear_fl = df_linear %>%
  ggplot(aes(var2, var1, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", na.value = "grey",
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name = NULL) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  geom_text(aes(var2, var1, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Pearson (Filtering)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1, 
                                   face = "italic"),
        axis.text.y = element_text(size = 12, face = "italic"),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 15),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") +
  coord_fixed()

heat_linear_fl

## -----------------------------------------------------------------------------
corr_dist = res_dist$dcorr_fl
cooccur_dist = res_dist$mat_cooccur

# Filter by co-occurrence
overlap = 10
corr_dist[cooccur_dist < overlap] = 0

df_dist = data.frame(get_upper_tri(corr_dist)) %>%
  rownames_to_column("var1") %>%
  pivot_longer(cols = -var1, names_to = "var2", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(value = round(value, 2))

tax_name = sort(union(df_dist$var1, df_dist$var2))
df_dist$var1 = factor(df_dist$var1, levels = tax_name)
df_dist$var2 = factor(df_dist$var2, levels = tax_name)

heat_dist_fl = df_dist %>%
  ggplot(aes(var2, var1, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", na.value = "grey",
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name = NULL) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  geom_text(aes(var2, var1, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Distance (Filtering)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1, 
                                   face = "italic"),
        axis.text.y = element_text(size = 12, face = "italic"),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 15),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") +
  coord_fixed()

heat_dist_fl

## -----------------------------------------------------------------------------
# Select subjects from "CE" and "NE"
tse1 = tse[, tse$region == "CE"]
tse2 = tse[, tse$region == "NE"]

# Rename samples to ensure there is an overlap of samples between CE and NE
colnames(tse1) = paste0("Sample-", seq_len(ncol(tse1)))
colnames(tse2) = paste0("Sample-", seq_len(ncol(tse2)))

print(tse1)
print(tse2)

## -----------------------------------------------------------------------------
set.seed(123)
# Linear relationships
res_linear = secom_linear(data = list(CE = tse1, NE = tse2), 
                          assay_name = c("counts", "counts"),
                          tax_level = c("Phylum", "Phylum"), pseudo = 0, 
                          prv_cut = 0.5, lib_cut = 1000, corr_cut = 0.5, 
                          wins_quant = c(0.05, 0.95), method = "pearson", 
                          soft = FALSE, thresh_len = 20, n_cv = 10, 
                          thresh_hard = 0.3, max_p = 0.005, n_cl = 2)

# Nonlinear relationships
res_dist = secom_dist(data = list(CE = tse1, NE = tse2),
                      assay_name = c("counts", "counts"),
                      tax_level = c("Phylum", "Phylum"), pseudo = 0, 
                      prv_cut = 0.5, lib_cut = 1000, corr_cut = 0.5, 
                      wins_quant = c(0.05, 0.95), R = 1000, 
                      thresh_hard = 0.3, max_p = 0.005, n_cl = 2)

## ----fig.width=8, fig.height=8------------------------------------------------
corr_linear = res_linear$corr_th
cooccur_linear = res_linear$mat_cooccur

# Filter by co-occurrence
overlap = 10
corr_linear[cooccur_linear < overlap] = 0

df_linear = data.frame(get_upper_tri(corr_linear)) %>%
  rownames_to_column("var1") %>%
  pivot_longer(cols = -var1, names_to = "var2", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(var2 = gsub("\\...", " - ", var2),
         value = round(value, 2))

tax_name = sort(union(df_linear$var1, df_linear$var2))
df_linear$var1 = factor(df_linear$var1, levels = tax_name)
df_linear$var2 = factor(df_linear$var2, levels = tax_name)
txt_color = ifelse(grepl("CE", tax_name), "#1B9E77", "#D95F02")

heat_linear_th = df_linear %>%
  ggplot(aes(var2, var1, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "grey", midpoint = 0, limit = c(-1,1), 
                       space = "Lab", name = NULL) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  geom_text(aes(var2, var1, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Pearson (Thresholding)") +
  theme_bw() +
  geom_vline(xintercept = 6.5, color = "blue", linetype = "dashed") +
  geom_hline(yintercept = 6.5, color = "blue", linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1, 
                                   face = "italic", color = txt_color),
        axis.text.y = element_text(size = 12, face = "italic", 
                                   color = txt_color),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 15),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") +
  coord_fixed()

heat_linear_th

## ----fig.width=8, fig.height=8------------------------------------------------
corr_linear = res_linear$corr_th
cooccur_linear = res_linear$mat_cooccur

# Filter by co-occurrence
overlap = 10
corr_linear[cooccur_linear < overlap] = 0

df_linear = data.frame(get_upper_tri(corr_linear)) %>%
  rownames_to_column("var1") %>%
  pivot_longer(cols = -var1, names_to = "var2", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(var2 = gsub("\\...", " - ", var2),
         value = round(value, 2))

tax_name = sort(union(df_linear$var1, df_linear$var2))
df_linear$var1 = factor(df_linear$var1, levels = tax_name)
df_linear$var2 = factor(df_linear$var2, levels = tax_name)
txt_color = ifelse(grepl("CE", tax_name), "#1B9E77", "#D95F02")

heat_linear_fl = df_linear %>%
  ggplot(aes(var2, var1, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "grey", midpoint = 0, limit = c(-1,1), 
                       space = "Lab", name = NULL) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  geom_text(aes(var2, var1, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Pearson (Filtering)") +
  theme_bw() +
  geom_vline(xintercept = 6.5, color = "blue", linetype = "dashed") +
  geom_hline(yintercept = 6.5, color = "blue", linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1, 
                                   face = "italic", color = txt_color),
        axis.text.y = element_text(size = 12, face = "italic", 
                                   color = txt_color),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 15),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") +
  coord_fixed()

heat_linear_fl

## ----fig.width=8, fig.height=8------------------------------------------------
corr_dist = res_dist$dcorr_fl
cooccur_dist = res_dist$mat_cooccur

# Filter by co-occurrence
overlap = 10
corr_dist[cooccur_dist < overlap] = 0

df_dist = data.frame(get_upper_tri(corr_dist)) %>%
  rownames_to_column("var1") %>%
  pivot_longer(cols = -var1, names_to = "var2", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(var2 = gsub("\\...", " - ", var2),
         value = round(value, 2))

tax_name = sort(union(df_dist$var1, df_dist$var2))
df_dist$var1 = factor(df_dist$var1, levels = tax_name)
df_dist$var2 = factor(df_dist$var2, levels = tax_name)
txt_color = ifelse(grepl("CE", tax_name), "#1B9E77", "#D95F02")

heat_dist_fl = df_dist %>%
  ggplot(aes(var2, var1, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "grey", midpoint = 0, limit = c(-1,1), 
                       space = "Lab", name = NULL) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  geom_text(aes(var2, var1, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Distance (Filtering)") +
  theme_bw() +
  geom_vline(xintercept = 6.5, color = "blue", linetype = "dashed") +
  geom_hline(yintercept = 6.5, color = "blue", linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1, 
                                   face = "italic", color = txt_color),
        axis.text.y = element_text(size = 12, face = "italic", 
                                   color = txt_color),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 15),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") +
  coord_fixed()

heat_dist_fl

## ----sessionInfo, message = FALSE, warning = FALSE, comment = NA--------------
sessionInfo()

