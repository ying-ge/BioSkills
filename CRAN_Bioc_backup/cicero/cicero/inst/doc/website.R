## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----getPackage, eval=FALSE---------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#   BiocManager::install("cicero")

## ----eval = FALSE-------------------------------------------------------------
#  BiocManager::install(cole-trapnell-lab/cicero)

## ----Load, message=FALSE------------------------------------------------------
 library(cicero)

## -----------------------------------------------------------------------------
data(cicero_data)

## ----eval=TRUE----------------------------------------------------------------
input_cds <- make_atac_cds(cicero_data, binarize = TRUE)

## ----eval=TRUE----------------------------------------------------------------
set.seed(2017)
input_cds <- detectGenes(input_cds)
input_cds <- estimateSizeFactors(input_cds)
input_cds <- reduceDimension(input_cds, max_components = 2, num_dim=6,
                        reduction_method = 'tSNE', norm_method = "none")

## ----eval=FALSE---------------------------------------------------------------
#  tsne_coords <- t(reducedDimA(input_cds))
#  row.names(tsne_coords) <- row.names(pData(input_cds))
#  cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = tsne_coords)

## ----eval=FALSE---------------------------------------------------------------
#  data("human.hg19.genome")
#  sample_genome <- subset(human.hg19.genome, V1 == "chr18")
#  sample_genome$V2[1] <- 10000000
#  conns <- run_cicero(cicero_cds, sample_genome, sample_num = 2) # Takes a few minutes to run
#  head(conns)

## ----fig.width = 7, fig.height = 4, fig.align='center', eval=FALSE------------
#  data(gene_annotation_sample)
#  plot_connections(conns, "chr18", 8575097, 8839855,
#                   gene_model = gene_annotation_sample,
#                   coaccess_cutoff = .25,
#                   connection_width = .5,
#                   collapseTranscripts = "longest" )

## ----eval=FALSE---------------------------------------------------------------
#  chia_conns <-  data.frame(Peak1 = c("chr18_10000_10200", "chr18_10000_10200",
#                                      "chr18_49500_49600"),
#                            Peak2 = c("chr18_10600_10700", "chr18_111700_111800",
#                                      "chr18_10600_10700"))
#  
#  conns$in_chia <- compare_connections(conns, chia_conns)

## ----eval=FALSE---------------------------------------------------------------
#  conns$in_chia_100 <- compare_connections(conns, chia_conns, maxgap=100)
#  
#  head(conns)

## ----eval=FALSE---------------------------------------------------------------
#  # Add a column of 1s called "coaccess"
#  chia_conns <-  data.frame(Peak1 = c("chr18_10000_10200", "chr18_10000_10200",
#                                      "chr18_49500_49600"),
#                            Peak2 = c("chr18_10600_10700", "chr18_111700_111800",
#                                      "chr18_10600_10700"),
#                            coaccess = c(1, 1, 1))
#  
#  plot_connections(conns, "chr18", 10000, 112367,
#                   gene_model = gene_annotation_sample,
#                   coaccess_cutoff = 0,
#                   connection_width = .5,
#                   comparison_track = chia_conns,
#                   include_axis_track = FALSE,
#                   collapseTranscripts = "longest")

## ----eval=FALSE---------------------------------------------------------------
#  CCAN_assigns <- generate_ccans(conns)
#  
#  head(CCAN_assigns)

## ----eval=FALSE---------------------------------------------------------------
#  
#  #### Add a column for the pData table indicating the gene if a peak is a promoter ####
#  # Create a gene annotation set that only marks the transcription start sites of
#  # the genes. We use this as a proxy for promoters.
#  # To do this we need the first exon of each transcript
#  pos <- subset(gene_annotation_sample, strand == "+")
#  pos <- pos[order(pos$start),]
#  pos <- pos[!duplicated(pos$transcript),] # remove all but the first exons per transcript
#  pos$end <- pos$start + 1 # make a 1 base pair marker of the TSS
#  
#  neg <- subset(gene_annotation_sample, strand == "-")
#  neg <- neg[order(neg$start, decreasing = TRUE),]
#  neg <- neg[!duplicated(neg$transcript),] # remove all but the first exons per transcript
#  neg$start <- neg$end - 1
#  
#  gene_annotation_sub <- rbind(pos, neg)
#  
#  # Make a subset of the TSS annotation columns containing just the coordinates
#  # and the gene name
#  gene_annotation_sub <- gene_annotation_sub[,c(1:3, 8)]
#  
#  # Rename the gene symbol column to "gene"
#  names(gene_annotation_sub)[4] <- "gene"
#  
#  input_cds <- annotate_cds_by_site(input_cds, gene_annotation_sub)
#  
#  head(fData(input_cds))
#  
#  
#  #### Generate gene activity scores ####
#  # generate unnormalized gene activity matrix
#  unnorm_ga <- build_gene_activity_matrix(input_cds, conns)
#  
#  # remove any rows/columns with all zeroes
#  unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, !Matrix::colSums(unnorm_ga) == 0]
#  
#  # make a list of num_genes_expressed
#  num_genes <- pData(input_cds)$num_genes_expressed
#  names(num_genes) <- row.names(pData(input_cds))
#  
#  # normalize
#  cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)
#  
#  # if you had two datasets to normalize, you would pass both:
#  # num_genes should then include all cells from both sets
#  unnorm_ga2 <- unnorm_ga
#  cicero_gene_activities <- normalize_gene_activities(list(unnorm_ga, unnorm_ga2), num_genes)
#  

## ----eval=TRUE----------------------------------------------------------------
data("cicero_data")
input_cds <- make_atac_cds(cicero_data)

# Add some cell meta-data
data("cell_data")
pData(input_cds) <- cbind(pData(input_cds), cell_data[row.names(pData(input_cds)),])
pData(input_cds)$cell <- NULL

agg_cds <- aggregate_nearby_peaks(input_cds, distance = 10000)
agg_cds <- detectGenes(agg_cds)
agg_cds <- estimateSizeFactors(agg_cds)
agg_cds <- estimateDispersions(agg_cds)

## ----eval=TRUE----------------------------------------------------------------
# This takes a few minutes to run
diff_timepoint <- differentialGeneTest(agg_cds,
                      fullModelFormulaStr="~timepoint + num_genes_expressed")

# We chose a very high q-value cutoff because there are
# so few sites in the sample dataset, in general a q-value
# cutoff in the range of 0.01 to 0.1 would be appropriate
ordering_sites <- row.names(subset(diff_timepoint, qval < .5))
length(ordering_sites)

## ----fig.show='hold', eval=FALSE----------------------------------------------
#  plot_pc_variance_explained(agg_cds, return_all = FALSE) #Choose 2 PCs
#  agg_cds <- reduceDimension(agg_cds,
#                                max_components = 2,
#                                norm_method = 'log',
#                                num_dim = 2,
#                                reduction_method = 'tSNE',
#                                verbose = TRUE)
#  
#  agg_cds <- clusterCells(agg_cds, verbose = FALSE)
#  
#  plot_cell_clusters(agg_cds, color_by = 'as.factor(Cluster)') + theme(text = element_text(size=8))
#  clustering_DA_sites <- differentialGeneTest(agg_cds[1:10,], #Subset for example only
#                                               fullModelFormulaStr = '~Cluster')
#  
#  # Not run because using Option 1 to continue
#  # ordering_sites <-
#  #  row.names(clustering_DA_sites)[order(clustering_DA_sites$qval)][1:1000]
#  

## ----eval=TRUE----------------------------------------------------------------
agg_cds <- setOrderingFilter(agg_cds, ordering_sites)

## ----fig.align='center', fig.height=4, fig.width=4, eval=TRUE-----------------
agg_cds <- reduceDimension(agg_cds, max_components = 2,
          residualModelFormulaStr="~num_genes_expressed",
          reduction_method = 'DDRTree')

# skipped due to monocle bug
# agg_cds <- suppressWarnings(orderCells(agg_cds))

plot_cell_trajectory(agg_cds, color_by = "timepoint")

## ----fig.align='center', fig.height=4, fig.width=4, eval=TRUE-----------------
# skipped due to monocle bug
#plot_cell_trajectory(agg_cds, color_by = "State")

## ----fig.align='center', fig.height=4, fig.width=4, eval=TRUE-----------------
# skipped due to monocle bug
#agg_cds <- suppressWarnings(orderCells(agg_cds, root_state = 4))
#plot_cell_trajectory(agg_cds, color_by = "Pseudotime")

## ----fig.align='center', fig.height=4, fig.width=4, eval=TRUE-----------------
# skipped due to monocle bug
#pData(input_cds)$Pseudotime <- pData(agg_cds)[colnames(input_cds),]$Pseudotime
#pData(input_cds)$State <- pData(agg_cds)[colnames(input_cds),]$State

## ----fig.width = 3, fig.height = 4, fig.align='center', eval=TRUE-------------
# skipped due to monocle bug
#input_cds_lin <- input_cds[,row.names(subset(pData(input_cds), State  != 5))]

#plot_accessibility_in_pseudotime(input_cds_lin[c("chr18_38156577_38158261", 
#                                                 "chr18_48373358_48374180", 
#                                                 "chr18_60457956_60459080")])

## ----eval=TRUE----------------------------------------------------------------
# skipped due to monocle bug
#pData(input_cds_lin)$cell_subtype <- cut(pData(input_cds_lin)$Pseudotime, 10)
#binned_input_lin <- aggregate_by_cell_bin(input_cds_lin, "cell_subtype")

## ----eval=TRUE----------------------------------------------------------------
# skipped due to monocle bug
#diff_test_res <- differentialGeneTest(binned_input_lin[1:10,], #Subset for example only
#    fullModelFormulaStr="~sm.ns(Pseudotime, df=3) + sm.ns(num_genes_expressed, df=3)",
#    reducedModelFormulaStr="~sm.ns(num_genes_expressed, df=3)", cores=1)

#head(diff_test_res)

## -----------------------------------------------------------------------------
citation("cicero")

## -----------------------------------------------------------------------------
sessionInfo()

