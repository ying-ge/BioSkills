## ----requirements, echo=FALSE-------------------------------------------------
if (!require(flowWorkspaceData)) {
  stop("Cannot build the vignettes without 'flowWorkspaceData'")
}

## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
opts_chunk$set(out.extra='style="display:block; margin: auto"', fig.align="center", message = FALSE, warning = FALSE)

## ----load-flowWorkspace, echo=F-----------------------------------------------
library(flowWorkspace)

## ----load-xml, eval=TRUE------------------------------------------------------
flowDataPath <- system.file("extdata", package = "flowWorkspaceData")
wsfile <- list.files(flowDataPath, pattern="manual.xml",full = TRUE)
wsfile

## ----open_flowjo_xml, eval=F--------------------------------------------------
#  library(CytoML)
#  ws <- open_flowjo_xml(wsfile)

## ----flowjo_to_gatingset, eval=F----------------------------------------------
#  gs <- flowjo_to_gatingset(ws, name= "T-cell", subset =1)

## ----load_gs_manual, echo = FALSE---------------------------------------------
gs <- load_gs(file.path(flowDataPath,"gs_manual"))

## ----plot-manual-GatingHierarchy----------------------------------------------
gh <- gs[[1]]
plot(gh)

## ----plot-manual-gates, fig.width = 9-----------------------------------------
library(ggcyto)
autoplot(gh)

## ----gatingTemplate, eval = T-------------------------------------------------
library(openCyto)
library(data.table)
gtFile <- system.file("extdata/gating_template/tcell.csv", package = "openCyto")
dtTemplate <- fread(gtFile)
dtTemplate

## ----gatingTemplate-nonDebris, eval = T---------------------------------------
dtTemplate[1,]

## ----gatingTemplate-singlets, eval = T----------------------------------------
dtTemplate[2,]

## ----gatingTemplate-lympth, eval = T------------------------------------------
dtTemplate[3,]

## ----gatingTemplate-cd3, eval = T---------------------------------------------
dtTemplate[4,]

## ----gatingTemplate-cd4cd8, eval = T------------------------------------------
dtTemplate[5,]

## ----gatingTemplate-expand, echo = F, results = F-----------------------------
expanded <- openCyto:::.preprocess_csv(dtTemplate)
rownames(expanded) <- NULL

## ----gatingTemplate-expand1, echo = F-----------------------------------------
expanded[5:6,]

## ----gatingTemplate-expand2, echo = F-----------------------------------------
expanded[7:10,]

## ----load-gt, eval = T--------------------------------------------------------
gt_tcell <- gatingTemplate(gtFile)
gt_tcell

## ----plot-gt, eval = T--------------------------------------------------------
plot(gt_tcell)

## ----load-fcs-----------------------------------------------------------------
fcsFiles <- list.files(pattern = "CytoTrol", flowDataPath, full = TRUE)
cs  <- load_cytoset_from_fcs(fcsFiles)
cf <- realize_view(cs[[1]])
gs <- GatingSet(cs)
gs

## ----compensate---------------------------------------------------------------
compMat <- gh_get_compensations(gh)
compensate(gs, compMat)

## ----compensate_plot, echo = F, fig.width = 4, fig.height = 4-----------------
sub_chnl <- c("V545-A","V450-A")
cf <- cf[,sub_chnl]
cf_comp <- realize_view(gh_pop_get_data(gs[[1]])[,sub_chnl])
cs <- cytoset(list(cf = cf, cf_comp = cf_comp))
#transform data to better visualize the compensation effect
transform(cs, estimateLogicle(cf,sub_chnl))
gridExtra::grid.arrange(as.ggplot(autoplot(cs[[1]], "V545", "V450")), as.ggplot(autoplot(cs[[2]], "V545", "V450")), nrow = 1)
# Need to grab the post-compensation, pre-transformation cytoframe before next steps
cf_comp <- realize_view(gh_pop_get_data(gs[[1]])[,sub_chnl])

## ----transformation, eval = T-------------------------------------------------
chnls <- parameters(compMat)
trans <- estimateLogicle(gs[[1]], channels = chnls)
gs <- transform(gs, trans)

## ----transformation_plot, echo = F, fig.width = 5, fig.height = 5-------------
cf_trans <- gh_pop_get_data(gs[[1]])[,sub_chnl[1]]
cf_comp <- cf_comp[,sub_chnl[1]]
p1 <- as.ggplot(autoplot(cf_comp, "V545"))
p2 <- as.ggplot(autoplot(cf_trans, "V545"))
plot(gridExtra::arrangeGrob(p1,p2))

## ----gating, eval = TRUE------------------------------------------------------
gt_gating(gt_tcell, gs)

## ----gating_par, eval = FALSE-------------------------------------------------
#  gt_gating(gt_tcell, gs, mc.cores=2, parallel_type = "multicore")

## ----plot_afterGating---------------------------------------------------------
plot(gs[[1]])

## ----hideGate, results = "hide"-----------------------------------------------
nodesToHide <- c("cd8+", "cd4+"
				, "cd4-cd8-", "cd4+cd8+"
				, "cd4+cd8-/HLA+", "cd4+cd8-/CD38+"
				, "cd4-cd8+/HLA+", "cd4-cd8+/CD38+"
				, "CD45_neg/CCR7_gate", "cd4+cd8-/CD45_neg"
				, "cd4-cd8+/CCR7+", "cd4-cd8+/CD45RA+"
				)
lapply(nodesToHide, function(thisNode) gs_pop_set_visibility(gs, thisNode, FALSE))

## ----rename, results = "hide"-------------------------------------------------
gs_pop_set_name(gs, "cd4+cd8-", "cd4")
gs_pop_set_name(gs, "cd4-cd8+", "cd8")

## ----plot_afterHiding---------------------------------------------------------
plot(gs[[1]])

## ----plotGate_autoGate, fig.width = 9-----------------------------------------
autoplot(gs[[1]])

## ----gt_add_gating_method-----------------------------------------------------
gs_add_gating_method(gs, alias = "non-activated cd4",
                         pop = "--",
                         parent = "cd4",
                         dims = "CD38,HLA",
                         gating_method = "tailgate")
plot(gs[[1]])

## ----gs_remove_gating_method--------------------------------------------------
gs_remove_gating_method(gs)
plot(gs[[1]])

