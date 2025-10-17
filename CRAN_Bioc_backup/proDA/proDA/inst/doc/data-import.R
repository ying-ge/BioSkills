## ----include = FALSE------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -------------------------------------------------------------------------------------------------
system.file("extdata/proteinGroups.txt", 
            package = "proDA", mustWork = TRUE)

## -------------------------------------------------------------------------------------------------
full_data <- read.delim(
    system.file("extdata/proteinGroups.txt", 
                package = "proDA", mustWork = TRUE),
    stringsAsFactors = FALSE
)

head(colnames(full_data))

## -------------------------------------------------------------------------------------------------
# I use a regular expression (regex) to select the intensity columns
intensity_colnames <- grep("^LFQ\\.intensity\\.", colnames(full_data), value=TRUE)

# Create matrix which only contains the intensity columns
data <- as.matrix(full_data[, intensity_colnames])
colnames(data) <- sub("^LFQ\\.intensity\\.", "", intensity_colnames)
# Code missing values explicitly as NA
data[data == 0] <- NA
# log transformation to account for mean-variance relation
data <- log2(data)
# Overview of data
data[1:7, 1:6]
# Set rownames after showing data, because they are so long
rownames(data) <- full_data$Protein.IDs

## -------------------------------------------------------------------------------------------------
annotation_df <- data.frame(
    Condition = sub("\\.\\d+", "", sub("^LFQ\\.intensity\\.", 
                                       "", intensity_colnames)),
    Replicate = as.numeric(sub("^LFQ\\.intensity\\.[[:alnum:]]+\\.", 
                               "", intensity_colnames)),
    stringsAsFactors = FALSE, row.names = colnames(data)
)
head(annotation_df)

## ----eval=FALSE, include=TRUE---------------------------------------------------------------------
#  # Not Run
#  library(proDA)
#  fit <- proDA(data, design= annotation_df$Condition, col_data = annotation_df)
#  test_diff(fit, contrast = CG1407 - S2R)
#  # End Not Run

## ----include=FALSE--------------------------------------------------------------------------------
library(SummarizedExperiment)
library(MSnbase)

## -------------------------------------------------------------------------------------------------
library(SummarizedExperiment)
se <- SummarizedExperiment(SimpleList(LFQ=data), colData=annotation_df)
rowData(se) <- full_data[, c("Only.identified.by.site", 
                             "Reverse", "Potential.contaminant")]
se

## -------------------------------------------------------------------------------------------------
library(MSnbase)

fData <- AnnotatedDataFrame(full_data[, c("Only.identified.by.site", 
                                 "Reverse", "Potential.contaminant")])
rownames(fData) <- rownames(data)
ms <- MSnSet(data, pData=AnnotatedDataFrame(annotation_df), fData=fData)
ms

## ----eval=FALSE, include=TRUE---------------------------------------------------------------------
#  # Not Run
#  library(proDA)
#  fit <- proDA(se, design = ~ Condition - 1)
#  test_diff(fit, contrast = ConditionCG1407 - ConditionS2R)
#  # End Not Run

## ----include=FALSE--------------------------------------------------------------------------------
library(dplyr)
library(stringr)
library(readr)
library(tidyr)
library(tibble)

## -------------------------------------------------------------------------------------------------
library(dplyr)
library(stringr)
library(readr)
library(tidyr)
library(tibble)
# Or short 
# library(tidyverse)

## -------------------------------------------------------------------------------------------------
# The read_tsv function works faster and more reliable than read.delim
# But it sometimes needs help to identify the right type for each column,
# because it looks only at the first 1,000 elements. 
# Here, I explicitly define the `Reverse` column as a character column
full_data <- read_tsv(
    system.file("extdata/proteinGroups.txt", 
                package = "proDA", mustWork = TRUE),
    col_types = cols(Reverse = col_character())
)

full_data

## -------------------------------------------------------------------------------------------------
# I explicitly call `dplyr::select()` because there is a naming conflict
# between the tidyverse and BioConductor packages for `select()` function
tidy_data <- full_data %>%
    dplyr::select(ProteinID=Protein.IDs, starts_with("LFQ.intensity.")) %>%
    gather(Sample, Intensity, starts_with("LFQ.intensity.")) %>%
    mutate(Condition = str_match(Sample, 
                 "LFQ\\.intensity\\.([[:alnum:]]+)\\.\\d+")[,2]) %>%
    mutate(Replicate = as.numeric(str_match(Sample, 
                 "LFQ\\.intensity\\.[[:alnum:]]+\\.(\\d+)")[,2])) %>%
    mutate(SampleName = paste0(Condition, ".", Replicate))

tidy_data

## -------------------------------------------------------------------------------------------------
data <- tidy_data %>%
    mutate(Intensity = ifelse(Intensity == 0, NA, log2(Intensity))) %>%
    dplyr::select(ProteinID, SampleName, Intensity) %>%
    spread(SampleName, Intensity) %>%
    column_to_rownames("ProteinID") %>%
    as.matrix()

data[1:4, 1:7]

annotation_df <- tidy_data %>%
    dplyr::select(SampleName, Condition, Replicate) %>%
    distinct() %>%
    arrange(Condition, Replicate) %>%
    as.data.frame() %>%
    column_to_rownames("SampleName")

annotation_df

## -------------------------------------------------------------------------------------------------
library(SummarizedExperiment)
se <- SummarizedExperiment(SimpleList(LFQ=data), colData=annotation_df)
rowData(se) <- full_data[, c("Only.identified.by.site", 
                             "Reverse", "Potential.contaminant")]
se

## -------------------------------------------------------------------------------------------------
library(MSnbase)

fData <- AnnotatedDataFrame(full_data[, c("Only.identified.by.site", 
                                 "Reverse", "Potential.contaminant")])
rownames(fData) <- rownames(data)
ms <- MSnSet(data, pData=AnnotatedDataFrame(annotation_df), fData=fData)
ms

## ----eval=FALSE, include=TRUE---------------------------------------------------------------------
#  # Not Run
#  library(proDA)
#  fit <- proDA(se, design = ~ Condition - 1)
#  test_diff(fit, contrast = ConditionCG1407 - ConditionS2R)
#  # End Not Run

## -------------------------------------------------------------------------------------------------
library(DEP)

full_data <- read.delim(
    system.file("extdata/proteinGroups.txt", 
                package = "proDA", mustWork = TRUE),
    stringsAsFactors = FALSE
)


exp_design <- data.frame(
   label =c("LFQ.intensity.CG1407.01", "LFQ.intensity.CG1407.02",  "LFQ.intensity.CG1407.03",  "LFQ.intensity.CG4676.01",  "LFQ.intensity.CG4676.02", "LFQ.intensity.CG4676.03", "LFQ.intensity.CG51963.01", "LFQ.intensity.CG51963.02", "LFQ.intensity.CG51963.03","LFQ.intensity.CG5620A.01", "LFQ.intensity.CG5620A.02", "LFQ.intensity.CG5620A.03", "LFQ.intensity.CG5620B.01","LFQ.intensity.CG5620B.02", "LFQ.intensity.CG5620B.03", "LFQ.intensity.CG5880.01", "LFQ.intensity.CG5880.02",  "LFQ.intensity.CG5880.03",  "LFQ.intensity.CG6017.01",  "LFQ.intensity.CG6017.02", "LFQ.intensity.CG6017.03", "LFQ.intensity.CG6618.01",  "LFQ.intensity.CG6618.02",  "LFQ.intensity.CG6618.03",  "LFQ.intensity.CG6627.01", "LFQ.intensity.CG6627.02", "LFQ.intensity.CG6627.03",  "LFQ.intensity.CG8314.01", "LFQ.intensity.CG8314.02",  "LFQ.intensity.CG8314.03", "LFQ.intensity.GsbPI.001", "LFQ.intensity.GsbPI.002",  "LFQ.intensity.GsbPI.003",  "LFQ.intensity.S2R.01", "LFQ.intensity.S2R.02", "LFQ.intensity.S2R.03"),
   condition = c("CG1407", "CG1407", "CG1407", "CG4676", "CG4676", "CG4676", "CG51963", "CG51963", "CG51963", "CG5620A", "CG5620A", "CG5620A", "CG5620B", "CG5620B", "CG5620B", "CG5880", "CG5880", "CG5880", "CG6017", "CG6017", "CG6017", "CG6618", "CG6618", "CG6618", "CG6627", "CG6627", "CG6627", "CG8314", "CG8314", "CG8314", "GsbPI", "GsbPI", "GsbPI", "S2R", "S2R", "S2R" ),
   replicate = rep(1:3, times=12),
   stringsAsFactors = FALSE
)

se <- import_MaxQuant(full_data, exp_design)
se
assay(se)[1:5, 1:5]

## ----eval=FALSE, include=TRUE---------------------------------------------------------------------
#  # Not Run
#  library(proDA)
#  fit <- proDA(se, design = ~ condition - 1)
#  # Here, we need to be specific, because DEP also has a test_diff method
#  proDA::test_diff(fit, contrast = conditionCG1407 - conditionS2R)
#  # End Not Run

## -------------------------------------------------------------------------------------------------
sessionInfo()

