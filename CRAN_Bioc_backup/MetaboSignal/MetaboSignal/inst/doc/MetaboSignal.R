## ----setup, echo=FALSE--------------------------------------------------------
knitr::opts_chunk$set(message=FALSE, fig.path='figures/')

## ----include = FALSE----------------------------------------------------------
library(MetaboSignal)

## ----message = FALSE, tidy = TRUE---------------------------------------------
MS_keggFinder(KEGG_database="organism", match = "rattus")

## ----tidy = TRUE, tidy.opts=list(indent = 4, width.cutoff = 50)---------------
MS_keggFinder(KEGG_database ="pathway", match = c("glycol", "inositol phosphate","insulin signal", "akt"), organism_code = "rno")

## ----tidy = TRUE--------------------------------------------------------------
metabo_paths <- c("rno00010","rno00562")
signaling_paths <- c("rno04910", "rno04151")

## ----tidy = TRUE, tidy.opts=list(indent = 4, width.cutoff = 50), results='asis',eval=FALSE----
#  MetaboSignal_table <- MS_keggNetwork(metabo_paths = metabo_paths,
#                                       signaling_paths = signaling_paths)
#  

## ----tidy = TRUE, tidy.opts=list(indent = 4, width.cutoff = 50)---------------
MetaboSignal_table <- MS_replaceNode(node1 = c("cpd:C00267", "cpd:C00221"), 
                                     node2 = "cpd:C00031", MetaboSignal_table)

## ----tidy = TRUE, tidy.opts=list(indent = 4, width.cutoff = 50), message = FALSE----
MS_findMappedNodes(nodes = c("cpd:C00267", "cpd:C00221", "cpd:C00031"),
                   MetaboSignal_table)

## ----tidy = TRUE, tidy.opts=list(indent = 4, width.cutoff = 50), message = FALSE----
## Get KEGG IDs
MS_convertGene(genes = c("303565", "65038", "309179"), organism_code = "rno", organism_name = "rat", output = "matrix")
MS_distances(MetaboSignal_table, organism_code = "rno", 
             source_genes = c("K01084", "K15909", "K11584"), 
             target_metabolites = "cpd:C00031", names = TRUE)

## ----tidy = TRUE, tidy.opts=list(indent = 4, width.cutoff = 50), eval = FALSE----
#  MS_shortestPathsNetwork(MetaboSignal_table, organism_code="rno",
#                          source_nodes = c("K01084", "K15909", "K11584"),
#                          target_nodes = "cpd:C00031", type = "bw",
#                          file_name = "MS")

