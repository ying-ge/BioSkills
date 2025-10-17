## ----setup, echo=FALSE--------------------------------------------------------
knitr::opts_chunk$set(message=FALSE, fig.path='figures/')

## ----message = FALSE, tidy = TRUE---------------------------------------------
## Load MetaboSignal
library(MetaboSignal)

## ----tidy = TRUE--------------------------------------------------------------
## Regulatory interactions
data("regulatory_interactions")
head(regulatory_interactions[, c(1,3,5)])

## KEGG metabolic pathways
data("kegg_pathways")
head(kegg_pathways[, -2])

## KEGG signaling pathways
tail(kegg_pathways[, -2])

## ----tidy = TRUE, eval=FALSE--------------------------------------------------
#  ## Get IDs of metabolic and signaling human pathways
#  hsa_paths <- MS_getPathIds(organism_code = "hsa")

## ----tidy = TRUE, tidy.opts=list(indent = 4, width.cutoff = 50)---------------
## Create metabo_paths and signaling_paths vectors
metabo_paths <- kegg_pathways[kegg_pathways[, "Path_type"] == "metabolic", "Path_id"]

signaling_paths <- kegg_pathways[kegg_pathways[, "Path_type"] == "signaling", "Path_id"]

## ----tidy = TRUE, tidy.opts=list(indent = 4, width.cutoff = 50), results='asis', eval=FALSE----
#  ## Build KEGG network (might take a while)
#  keggNet_example <- MS_keggNetwork(metabo_paths, signaling_paths, expand_genes = TRUE,
#                                    convert_entrez = TRUE)

## ----tidy = TRUE--------------------------------------------------------------
## See network format
head(keggNet_example)

## ----tidy = TRUE, tidy.opts=list(indent = 4, width.cutoff = 55)---------------
## Get all types of interaction
all_types <- unique(unlist(strsplit(keggNet_example[, "interaction_type"], "/")))
all_types <- gsub("k_", "", all_types)

## Select wanted interactions
wanted_types <- setdiff(all_types, c("unknown", "indirect-compound", "indirect-effect", 
                                     "dissociation", "state-change", "binding",
                                     "association"))
print(wanted_types) # interactions that will be retained

## Filter keggNet_example to retain only wanted interactions
wanted_types <- paste(wanted_types, collapse = "|")
keggNet_clean <- keggNet_example[grep(wanted_types, keggNet_example[, 3]), ]

## ----tidy = TRUE--------------------------------------------------------------
## Build regulatory network of TRRUST interactions
trrustNet_example <- MS2_ppiNetwork(datasets = "trrust")

## Build regulatory network of OmniPath interactions
omnipathNet_example <- MS2_ppiNetwork(datasets = "omnipath")

## Build regulatory network by merging OmniPath and TRRUST interactions
ppiNet_example <- MS2_ppiNetwork(datasets = "all")

## See network format
head(ppiNet_example)

## ----tidy = TRUE, tidy.opts=list(indent = 4, width.cutoff = 60), results='asis', eval=FALSE----
#  ## Merge networks
#  mergedNet_example <- MS2_mergeNetworks(keggNet_clean, ppiNet_example)

## ----message = FALSE, tidy = TRUE---------------------------------------------
## See network format
head(mergedNet_example)

