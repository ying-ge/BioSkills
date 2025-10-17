## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  out.width = "100%",
  error = TRUE
)
knitr::opts_knit$set(global.par = TRUE)

## ----include=FALSE--------------------------------------------------------------------------------
set.seed(1)
options(width = 100)
par(cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

## ----eval=FALSE-----------------------------------------------------------------------------------
#  if(!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("proDA")

## ----eval=FALSE-----------------------------------------------------------------------------------
#  # install.packages("devtools")
#  devtools::install_github("const-ae/proDA")

## ----quickstart-----------------------------------------------------------------------------------
# Load the package
library(proDA)
# Generate some dataset with known structure
syn_dataset <- generate_synthetic_data(n_proteins = 100, n_conditions = 2)

# The abundance matrix
syn_dataset$Y[1:5, ]

# Assignment of the samples to the two conditions
syn_dataset$groups

# Fit the probabilistic dropout model
fit <- proDA(syn_dataset$Y, design = syn_dataset$groups)

# Identify which proteins differ between Condition 1 and 2
test_diff(fit, `Condition_1` - `Condition_2`, sort_by = "pval", n_max = 5)

## -------------------------------------------------------------------------------------------------
system.file("extdata/proteinGroups.txt", package = "proDA", mustWork = TRUE)

## -------------------------------------------------------------------------------------------------
# Load the table into memory
maxquant_protein_table <- read.delim(
    system.file("extdata/proteinGroups.txt", package = "proDA", mustWork = TRUE),
    stringsAsFactors = FALSE
)

## -------------------------------------------------------------------------------------------------
# I use a regular expression (regex) to select the intensity columns
intensity_colnames <- grep("^LFQ\\.intensity\\.", colnames(maxquant_protein_table), value=TRUE)
head(intensity_colnames)


# Create the intensity matrix
abundance_matrix <- as.matrix(maxquant_protein_table[, intensity_colnames])
# Adapt column and row maxquant_protein_table
colnames(abundance_matrix) <- sub("^LFQ\\.intensity\\.", "", intensity_colnames)
rownames(abundance_matrix) <- maxquant_protein_table$Protein.IDs
# Print some rows of the matrix with short names so they fit on the screen
abundance_matrix[46:48, 1:6]

## -------------------------------------------------------------------------------------------------
abundance_matrix[abundance_matrix == 0] <- NA
abundance_matrix <- log2(abundance_matrix)
abundance_matrix[46:48, 1:6]


## ----qc-mis_barplot,  out.width="80%", fig.height=4, fig.align="center"---------------------------

barplot(colSums(is.na(abundance_matrix)),
        ylab = "# missing values",
        xlab = "Sample 1 to 36")

## ----qc-raw_boxplot, out.width="80%", fig.height=4, fig.align="center"----------------------------
boxplot(abundance_matrix,
        ylab = "Intensity Distribution",
        xlab = "Sample 1 to 36")

## -------------------------------------------------------------------------------------------------
normalized_abundance_matrix <- median_normalization(abundance_matrix)

## -------------------------------------------------------------------------------------------------
da <- dist_approx(normalized_abundance_matrix)

## ----sample_dist, out.width="60%", fig.align="center"---------------------------------------------
# This chunk only works if pheatmap is installed
# install.packages("pheatmap")
sel <- c(1:3,  # CG1407
         7:9,  # CG59163
         22:24)# CG6618

plot_mat <- as.matrix(da$mean)[sel, sel]
# Remove diagonal elements, so that the colorscale is not distorted
plot_mat[diag(9) == 1] <- NA
# 95% conf interval is approx `sd * 1.96`
uncertainty <- matrix(paste0(" ± ",round(as.matrix(da$sd * 1.96)[sel, sel], 1)), nrow=9)
pheatmap::pheatmap(plot_mat, 
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   display_numbers= uncertainty,
                   number_color = "black")

## -------------------------------------------------------------------------------------------------
# The best way to create this data.frame depends on the column naming scheme
sample_info_df <- data.frame(name = colnames(normalized_abundance_matrix),
                             stringsAsFactors = FALSE)
sample_info_df$condition <- substr(sample_info_df$name, 1, nchar(sample_info_df$name)  - 3)
sample_info_df$replicate <- as.numeric(
  substr(sample_info_df$name, nchar(sample_info_df$name)  - 1, 20)
)
sample_info_df

## -------------------------------------------------------------------------------------------------
fit <- proDA(normalized_abundance_matrix, design = ~ condition, 
             col_data = sample_info_df, reference_level = "S2R")
fit

## -------------------------------------------------------------------------------------------------
# Equivalent to feature_parameters(fit)
fit$feature_parameters

## ----protein_dist, out.width="60%", fig.align="center"--------------------------------------------
# This chunk only works if pheatmap is installed
# install.packages("pheatmap")
pheatmap::pheatmap(dist_approx(fit[1:20, 1:3], by_sample = FALSE)$mean)

## -------------------------------------------------------------------------------------------------
# Test which proteins differ between condition CG1407 and S2R
# S2R is the default contrast, because it was specified as the `reference_level`
test_res <- test_diff(fit, "conditionCG1407")
test_res

## -------------------------------------------------------------------------------------------------
sessionInfo()

