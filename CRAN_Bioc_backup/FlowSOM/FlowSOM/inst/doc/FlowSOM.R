### R code from vignette source 'FlowSOM.Rnw'

###################################################
### code chunk number 1: style-Sweave
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: FlowSOM.Rnw:54-56
###################################################
# install.packages("BiocManager")
# BiocManager::install("FlowSOM")


###################################################
### code chunk number 3: FlowSOM.Rnw:68-78
###################################################
fileName <- system.file("extdata", "68983.fcs", package="FlowSOM")
ff <- flowCore::read.FCS(fileName)

# Compensation
comp <- flowCore::keyword(ff)[["SPILL"]]
ff <- flowWorkspace::compensate(ff, comp)

# Transformation
transformList <- flowCore::estimateLogicle(ff, channels = colnames(comp))
ff <- flowWorkspace::transform(ff, transformList)


###################################################
### code chunk number 4: FlowSOM.Rnw:89-101
###################################################
set.seed(42)
library(FlowSOM)

fSOM <- FlowSOM(ff,
                # Input options:
                compensate = FALSE, 
                transform = FALSE,
                scale = FALSE,
                # SOM options:
                colsToUse = c(9, 12, 14:18), xdim = 7, ydim = 7,
                # Metaclustering options:
                nClus = 10)


###################################################
### code chunk number 5: FlowSOM.Rnw:104-106
###################################################
p <- PlotStars(fSOM, backgroundValues = fSOM$metaclustering)
print(p, newpage = FALSE)


###################################################
### code chunk number 6: FlowSOM.Rnw:112-114
###################################################
FlowSOMmary(fsom = fSOM,
            plotFile = "FlowSOMmary.pdf")


###################################################
### code chunk number 7: FlowSOM.Rnw:119-121
###################################################
head(GetClusters(fSOM))
head(GetMetaclusters(fSOM))


###################################################
### code chunk number 8: FlowSOM.Rnw:153-164
###################################################
set.seed(42)
library(flowCore)
library(FlowSOM)

fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
fSOM <- ReadInput(fileName, compensate = TRUE, transform = TRUE, 
                    toTransform = c(8:18), scale = TRUE)

ff <- suppressWarnings(flowCore::read.FCS(fileName))
fSOM <- ReadInput(ff, compensate = TRUE, transform = TRUE,
                    toTransform = c(8:18), scale = TRUE)


###################################################
### code chunk number 9: FlowSOM.Rnw:175-176
###################################################
str(fSOM, max.level = 2)


###################################################
### code chunk number 10: FlowSOM.Rnw:194-196
###################################################
fSOM <- BuildSOM(fSOM, colsToUse = c(9, 12, 14:18))
str(fSOM$map, max.level = 2)


###################################################
### code chunk number 11: FlowSOM.Rnw:206-208
###################################################
fSOM <- BuildMST(fSOM)
str(fSOM$MST, max.level =  1)


###################################################
### code chunk number 12: FlowSOM.Rnw:216-217
###################################################
PlotStars(fSOM)


###################################################
### code chunk number 13: FlowSOM.Rnw:220-221
###################################################
PlotStars(fSOM, view = "grid")


###################################################
### code chunk number 14: FlowSOM.Rnw:224-227
###################################################
set.seed(1)
tsne <- Rtsne::Rtsne(fSOM$map$codes, perplexity = 6)
PlotStars(fSOM, view = tsne$Y, maxNodeSize = 2)


###################################################
### code chunk number 15: FlowSOM.Rnw:233-235
###################################################
p <- PlotStars(fSOM, equalNodeSize = TRUE)
print(p, newpage = FALSE)


###################################################
### code chunk number 16: FlowSOM.Rnw:242-253
###################################################
wsp_file <- system.file("extdata", "gating.wsp", package = "FlowSOM")

# Specify the cell types of interest for assigning one label per cell
cell_types <- c("B cells",
                "gd T cells", "CD4 T cells", "CD8 T cells",
                "NK cells","NK T cells")

# Parse the FlowJo workspace
gatingResult <- GetFlowJoLabels(fileName, 
                                wsp_file,
                                cell_types = cell_types)


###################################################
### code chunk number 17: FlowSOM.Rnw:256-257
###################################################
PlotPies(fSOM, cellTypes = gatingResult$manual)


###################################################
### code chunk number 18: FlowSOM.Rnw:263-265
###################################################
p <- PlotMarker(fSOM, "PE-Cy5-A")
print(p, newpage = FALSE)


###################################################
### code chunk number 19: FlowSOM.Rnw:267-269
###################################################
p <- PlotMarker(fSOM, "CD19")
print(p, newpage = FALSE)


###################################################
### code chunk number 20: FlowSOM.Rnw:273-274
###################################################
PlotNumbers(fSOM)


###################################################
### code chunk number 21: FlowSOM.Rnw:278-283
###################################################
plot <- Plot2DScatters(fSOM, 
                       channelpairs = list(c("PE-Texas Red-A", "Pacific Blue-A")),
                       clusters = list(c(81, 82, 91, 92, 93)),
                       plotFile = NULL)
print(plot[[1]])


###################################################
### code chunk number 22: FlowSOM.Rnw:296-297
###################################################
metaClustering <- as.character(metaClustering_consensus(fSOM$map$codes,k = 7))


###################################################
### code chunk number 23: FlowSOM.Rnw:300-303
###################################################
PlotPies(fSOM,
         cellTypes=gatingResult$manual,
         backgroundValues = metaClustering)


###################################################
### code chunk number 24: FlowSOM.Rnw:306-307
###################################################
PlotLabels(fSOM, labels = metaClustering)


###################################################
### code chunk number 25: FlowSOM.Rnw:311-313
###################################################
metaClustering_perCell <- GetMetaclusters(fSOM, metaClustering)
table(metaClustering_perCell)


###################################################
### code chunk number 26: FlowSOM.Rnw:323-332
###################################################
# Look for CD8+ ab T cells
query <- c("PE-Cy7-A" = "high", #CD3
            "APC-Cy7-A" = "high", #TCRb
            "Pacific Blue-A" = "high") #CD8
query_res <- QueryStarPlot(fSOM, query, equalNodeSize = TRUE, plot = FALSE)

cellTypes <- factor(rep("Unlabeled", fSOM$map$nNodes), 
                    levels=c("Unlabeled", "CD8 T cells"))
cellTypes[query_res$selected] <- "CD8 T cells"


###################################################
### code chunk number 27: FlowSOM.Rnw:335-338
###################################################
p <- PlotStars(fSOM, backgroundValues=cellTypes, 
               backgroundColor=c("#FFFFFF00","#0000FF"))
print(p, newpage = FALSE)


###################################################
### code chunk number 28: FlowSOM.Rnw:344-391
###################################################
library(FlowSOM)
set.seed(1)

# Build FlowSom result
 fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
 ff <- flowCore::read.FCS(fileName)
 ff <- flowCore::compensate(ff, ff@description$SPILL)
 ff <- flowCore::transform(ff,
           flowCore::transformList(colnames(ff@description$SPILL),
                                     flowCore::logicleTransform()))
 flowSOM.res <- FlowSOM(ff, scale = TRUE, colsToUse = c(9, 12, 14:18),
                          nClus = 10)
   
# Create new data
# To illustrate the output, we here generate new fcs files (with more 
# cells in metaclusters 1 and 9).
# In practice you would not generate any new file but use your different
# files from your different groups
 for(i in 1:5){
   flowCore::write.FCS(ff[sample(nrow(ff), 1000), ], 
                       file = paste0("ff_tmp", i, ".fcs"))  
 }
 for(i in 6:10){
    flowCore::write.FCS(ff[c(sample(nrow(ff),500),
                          sample(which(GetMetaclusters(flowSOM.res) == 1), 250),
                          sample(which(GetMetaclusters(flowSOM.res) == 9), 250)), ], 
                     file = paste0("ff_tmp", i, ".fcs"))
 }                  
 
# Get the count matrix
 percentages <- GetFeatures(fsom = flowSOM.res, 
                            files = paste0("ff_tmp",1:10,".fcs"), 
                            type = "percentages")
  
# Perform the statistics
groups <- list("Group 1" = paste0("ff_tmp", 1:5, ".fcs"), 
               "Group 2" = paste0("ff_tmp", 6:10, ".fcs"))
MC_stats <- GroupStats(percentages[["metacluster_percentages"]], groups)
C_stats <- GroupStats(percentages[["cluster_percentages"]], groups)

# Process the fold changes vector
fold_changes <- C_stats["fold changes", ]
fold_changes <- factor(ifelse(fold_changes < -3, "Underrepresented compared to Group 1",
                       ifelse(fold_changes > 3, "Overrepresented compared to Group 1",
                       "--")), levels = c("--", "Underrepresented compared to Group 1",
                       "Overrepresented compared to Group 1"))
fold_changes[is.na(fold_changes)] <- "--"


###################################################
### code chunk number 29: FlowSOM.Rnw:394-413
###################################################
# Show in figure
## Fold change
gr_1 <- PlotStars(flowSOM.res, title = "Group 1", 
                  nodeSizes = C_stats["medians Group 1", ], 
                  refNodeSize = max(C_stats[c("medians Group 1", "medians Group 2"),]),
                  backgroundValues = fold_changes,
                  backgroundColors = c("white", "red", "blue"), 
                  list_insteadof_ggarrange = TRUE)
gr_2 <- PlotStars(flowSOM.res, title = "Group 2", 
                  nodeSizes = C_stats["medians Group 2", ], 
                  refNodeSize = max(C_stats[c("medians Group 1", "medians Group 2"),]),
                  backgroundValues = fold_changes,
                  backgroundColors = c("white", "red", "blue"), 
                  list_insteadof_ggarrange = TRUE)
p <- ggpubr::ggarrange(plotlist = list(gr_1$tree, gr_2$tree, 
                                       gr_2$starLegend, gr_2$backgroundLegend), 
                       ncol = 2, nrow = 2,
                       heights = c(4, 1))
print(p, newpage = FALSE)


###################################################
### code chunk number 30: FlowSOM.Rnw:416-422
###################################################
## p values
p <- PlotVariable(flowSOM.res, 
                  title = "Fold change group 1 vs. group 2",
                  variable = C_stats["fold changes", ])

print(p, newpage = FALSE)


###################################################
### code chunk number 31: FlowSOM.Rnw:425-449
###################################################
## volcano plot
p <- ggplot2::ggplot(data.frame("-log10 p values" = c(C_stats[4, ], 
                                                      MC_stats[4, ]), 
                                "log10 fold changes" = c(C_stats[7, ], 
                                                         MC_stats[7, ]), 
                                "feature" = c(colnames(C_stats), colnames(MC_stats)),
                                "metacluster" = c(as.character(flowSOM.res$metaclustering),
                                                  levels(flowSOM.res$metaclustering)),
                                "type" = c(rep("C", ncol(C_stats)),
                                           rep("MC", ncol(MC_stats))),
                                check.names = FALSE), 
                     ggplot2::aes(x = `log10 fold changes`, 
                                  y = `-log10 p values`,
                                  size = type,
                                  col = metacluster)) +
  ggplot2::xlim(-3, 3) +
  ggplot2::ylim(0, 3) +
  ggplot2::geom_point() +
  ggplot2::theme_minimal() +
  ggplot2::geom_vline(xintercept = 0, col = "darkgrey") +
  ggplot2::geom_hline(yintercept = -log10(0.05), col = "darkgrey") +
  ggplot2::scale_size_manual(values =c("C" = 1, "MC" = 2)) 

print(p, newpage = FALSE)


###################################################
### code chunk number 32: FlowSOM.Rnw:461-462
###################################################
sessionInfo()


