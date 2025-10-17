## ----include=FALSE, echo=FALSE------------------------------------------------
# date: "`r doc_date()`"
# "`r pkg_ver('BiocStyle')`"
# <style>
#     pre {
#     white-space: pre !important;
#     overflow-y: scroll !important;
#     height: 50vh !important;
#     }
# </style>

## ----annotate, echo=FALSE-------------------------------------------------------------------------
library(knitr)
options(width=100)
opts_chunk$set(echo       = TRUE,
               message    = TRUE,
               warning    = TRUE,
               eval       = TRUE)

## ----multimir_dbInfoVersions----------------------------------------------------------------------
library(multiMiR)
db.ver = multimir_dbInfoVersions()
db.ver

## ----multimir_switchDBVersion, echo=TRUE----------------------------------------------------------
vers_table <- multimir_dbInfoVersions()
vers_table

multimir_switchDBVersion(db_version = "2.0.0")

curr_vers  <- vers_table[1, "VERSION"]  # current version
multimir_switchDBVersion(db_version = curr_vers)

## ----multimir_dbTables----------------------------------------------------------------------------
db.tables = multimir_dbTables()
db.tables

## ----multimir_dbInfo------------------------------------------------------------------------------
db.info = multimir_dbInfo()
db.info

## ----multimir-tabletype---------------------------------------------------------------------------
predicted_tables()
validated_tables()
diseasedrug_tables()
reverse_table_lookup("targetscan")

## ----multimir_dbCount-----------------------------------------------------------------------------
db.count = multimir_dbCount()
db.count
apply(db.count[,-1], 2, sum)

## ----list_multimir--------------------------------------------------------------------------------
miRNAs   = list_multimir("mirna", limit = 10)
genes    = list_multimir("gene", limit = 10)
drugs    = list_multimir("drug", limit = 10)
diseases = list_multimir("disease", limit = 10)
# executes 2 separate queries, giving 20 results
head(miRNAs)
head(genes)
head(drugs)
head(diseases)

## ----biocworkflow, eval=TRUE----------------------------------------------------------------------
library(edgeR)
library(multiMiR)

# Load data
counts_file  <- system.file("extdata", "counts_table.Rds", package = "multiMiR")
strains_file <- system.file("extdata", "strains_factor.Rds", package = "multiMiR")
counts_table   <- readRDS(counts_file)
strains_factor <- readRDS(strains_file)

# Standard edgeR differential expression analysis
design <- model.matrix(~ strains_factor)

# Using trended dispersions
dge <- DGEList(counts = counts_table)
dge <- calcNormFactors(dge)
dge$samples$strains <- strains_factor
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

# Fit GLM model for strain effect
fit <- glmFit(dge, design)
lrt <- glmLRT(fit)

# Table of unadjusted p-values (PValue) and FDR values
p_val_DE_edgeR <- topTags(lrt, adjust.method = 'BH', n = Inf)

# Getting top differentially expressed miRNA's
top_miRNAs <- rownames(p_val_DE_edgeR$table)[1:10]

# Plug miRNA's into multiMiR and getting validated targets
multimir_results <- get_multimir(org     = 'mmu',
                                 mirna   = top_miRNAs,
                                 table   = 'validated',
                                 summary = TRUE)
head(multimir_results@data)

## ----Example1-------------------------------------------------------------------------------------
# The default is to search validated interactions in human
example1 <- get_multimir(mirna = 'hsa-miR-18a-3p', summary = TRUE)
names(example1)
# Check which types of associations were returned
table(example1@data$type)
# Detailed information of the validated miRNA-target interaction
head(example1@data)
# Which interactions are supported by Luciferase assay?
example1@data[grep("Luciferase", example1@data[, "experiment"]), ]
example1@summary[example1@summary[,"target_symbol"] == "KRAS",]

## ----Example2-------------------------------------------------------------------------------------
example2 <- get_multimir(disease.drug = 'cisplatin', table = 'disease.drug')
names(example2)
nrow(example2@data)
table(example2@data$type)
head(example2@data)

## ----Example3_part1-------------------------------------------------------------------------------
example3 <- get_multimir(org     = "mmu",
                         target  = "Gnb1",
                         table   = "predicted",
                         summary = TRUE,
                         predicted.cutoff      = 35,
                         predicted.cutoff.type = "p",
                         predicted.site        = "all")
names(example3)
table(example3@data$type)
head(example3@data)
head(example3@summary)

## ----Example3_part2-------------------------------------------------------------------------------
apply(example3@summary[, 6:13], 2, function(x) sum(x > 0))

## ----Example4_part1-------------------------------------------------------------------------------
example4 <- get_multimir(org     = 'hsa',
                         target  = c('AKT2', 'CERS6', 'S1PR3', 'SULF2'),
                         table   = 'predicted',
                         summary = TRUE,
                         predicted.cutoff.type = 'n',
                         predicted.cutoff      = 500000)

## ----Example4_part2-------------------------------------------------------------------------------
example4.counts <- addmargins(table(example4@summary[, 2:3]))
example4.counts <- example4.counts[-nrow(example4.counts), ]
example4.counts <- example4.counts[order(example4.counts[, 5], decreasing = TRUE), ]
head(example4.counts)

## ----Example5, echo=TRUE--------------------------------------------------------------------------
load(url("http://multimir.org/bladder.rda"))

## ----Example5_part2, eval=TRUE, echo=TRUE---------------------------------------------------------
# search all tables & top 10% predictions
example5 <- get_multimir(org     = "hsa",
                         mirna   = DE.miRNA.up,
                         target  = DE.entrez.dn,
                         table   = "all",
                         summary = TRUE,
                         predicted.cutoff.type = "p",
                         predicted.cutoff      = 10,
                         use.tibble = TRUE)

## ----eval = FALSE, include=FALSE------------------------------------------------------------------
#  ## Searching diana_microt ...
#  ## Searching elmmo ...
#  ## Searching microcosm ...
#  ## Searching miranda ...
#  ## Searching mirdb ...
#  ## Searching pictar ...
#  ## Searching pita ...
#  ## Searching targetscan ...
#  ## Searching mirecords ...
#  ## Searching mirtarbase ...
#  ## Searching tarbase ...
#  ## Searching mir2disease ...
#  ## Searching pharmaco_mir ...
#  ## Searching phenomir ...

## ----Example5_part3, eval=TRUE, echo=TRUE---------------------------------------------------------
table(example5@data$type)
result <- select(example5, keytype = "type", keys = "validated", columns = columns(example5))
unique_pairs <- 
    result[!duplicated(result[, c("mature_mirna_id", "target_entrez")]), ]

result

## ----eval=FALSE, include=FALSE--------------------------------------------------------------------
#  ##     database mature_mirna_acc mature_mirna_id target_symbol target_entrez
#  ## 1 mirtarbase     MIMAT0000087  hsa-miR-30a-5p          FDX1          2230
#  ## 2 mirtarbase     MIMAT0000087  hsa-miR-30a-5p        LIMCH1         22998
#  ## 3    tarbase     MIMAT0000087  hsa-miR-30a-5p          FDX1          2230
#  ## 4    tarbase     MIMAT0000424     hsa-miR-128          NEK2          4751
#  ## 5    tarbase     MIMAT0000087  hsa-miR-30a-5p        LIMCH1         22998
#  ##    target_ensembl               experiment          support_type pubmed_id
#  ## 1 ENSG00000137714               Proteomics Functional MTI (Weak)  18668040
#  ## 2 ENSG00000064042 pSILAC//Proteomics;Other Functional MTI (Weak)  18668040
#  ## 3 ENSG00000137714               Proteomics              positive
#  ## 4 ENSG00000117650               Microarray              positive
#  ## 5 ENSG00000064042               Proteomics              positive

## ----Example5_part4, eval=TRUE, echo=TRUE---------------------------------------------------------
mykeytype <- "disease_drug"

mykeys <- keys(example5, keytype = mykeytype)
mykeys <- mykeys[grep("bladder", mykeys, ignore.case = TRUE)]

result <- select(example5, keytype = "disease_drug", keys = mykeys,
                 columns = columns(example5))
result

## ----Example5_part4_fortext, echo=FALSE, include=FALSE, eval=TRUE---------------------------------
unique_pairs <- 
    result[!duplicated(apply(result[, c("mature_mirna_id", "disease_drug")], 2,
                             tolower)), ]

## ----eval=FALSE, include=FALSE--------------------------------------------------------------------
#  ##        database mature_mirna_acc mature_mirna_id target_symbol target_entrez
#  ## 18  mir2disease     MIMAT0000418  hsa-miR-23b-3p            NA            NA
#  ## 711    phenomir     MIMAT0000418  hsa-miR-23b-3p            NA            NA
#  ## 311    phenomir     MIMAT0000449 hsa-miR-146a-5p            NA            NA
#  ##     target_ensembl   disease_drug
#  ## 18              NA bladder cancer
#  ## 711             NA Bladder cancer
#  ## 311             NA Bladder cancer
#  ##                                               paper_pubmedID
#  ## 18  2007. Micro-RNA profiling in kidney and bladder cancers.
#  ## 711                                                 17826655
#  ## 311                                                 19127597

## ----Example5_part5, eval=TRUE, echo=TRUE---------------------------------------------------------
predicted <- select(example5, keytype = "type", keys = "predicted", 
                    columns = columns(example5))
length(unique(predicted$mature_mirna_id))
length(unique(predicted$target_entrez))
unique.pairs <- 
    unique(data.frame(miRNA.ID = as.character(predicted$mature_mirna_id),
                      target.Entrez = as.character(predicted$target_entrez)))
nrow(unique.pairs)

## ----Example5_part8, eval=TRUE, echo=TRUE---------------------------------------------------------
head(unique.pairs)

## ----Example5_part9, eval=TRUE, echo=TRUE---------------------------------------------------------
example5.split <- split(predicted, predicted$database)

## ----annodbi--------------------------------------------------------------------------------------
# On example4's result
columns(example4)
head(keys(example4))
keytypes(example4)
mykeys <- keys(example4)[1:4]
head(select(example4, keys = mykeys, 
            columns = c("database", "target_entrez")))

# On example3's result
columns(example3)
head(keys(example3))
keytypes(example3)
mykeys <- keys(example3)[1:4]
head(select(example3, keys = mykeys, 
            columns = c("database", "target_entrez", "score")))

# Search by gene on example4's result
columns(example4)
keytypes(example4)
head(keys(example4, keytype = "target_entrez"))
mykeys <- keys(example4, keytype = "target_entrez")[1]
head(select(example4, keys = mykeys, keytype = "target_entrez",
            columns = c("database", "target_entrez", "score")))

## ----Direct_query1--------------------------------------------------------------------------------
direct2 <- search_multimir(query = "describe diana_microt")
direct2

## ----Direct_query2--------------------------------------------------------------------------------
qry <- "SELECT m.mature_mirna_acc, m.mature_mirna_id, t.target_symbol,
               t.target_entrez, t.target_ensembl, i.experiment, i.support_type,
               i.pubmed_id 
        FROM mirna AS m INNER JOIN mirecords AS i INNER JOIN target AS t 
        ON (m.mature_mirna_uid=i.mature_mirna_uid and 
            i.target_uid=t.target_uid) 
        WHERE m.mature_mirna_id='hsa-miR-18a-3p'"
direct3 <- search_multimir(query = qry)
direct3

## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()
warnings()

