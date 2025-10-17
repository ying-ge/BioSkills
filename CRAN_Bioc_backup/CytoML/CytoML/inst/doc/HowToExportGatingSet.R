## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, results = "markup", message = FALSE, warning = FALSE)

## -----------------------------------------------------------------------------
library(ncdfFlow)
library(flowWorkspace)
library(CytoML)
dataDir <- system.file("extdata",package="flowWorkspaceData")
#load raw FCS
fs <- load_cytoset_from_fcs(file.path(dataDir,"CytoTrol_CytoTrol_1.fcs"))
gs <- GatingSet(fs)

## -----------------------------------------------------------------------------
#compensate
comp <- spillover(fs[[1]])[["SPILL"]]
chnls <- colnames(comp)
comp <- compensation(comp)
gs <- compensate(gs, comp)

#transform
trans <- flowjo_biexp_trans()
trans <- transformerList(chnls, trans)
gs <- transform(gs, trans)

## ----warning=FALSE------------------------------------------------------------
library(openCyto)
#load the original template for tcell panel
tbl <- data.table::fread(system.file("extdata/gating_template/tcell.csv", package = "openCyto"))
#modify some paramters to fit the current data range
tbl[5, gating_args:= "gate_range = c(1e3, 3e3)"]
tbl[c(8,11), gating_args:= "gate_range = c(2e3, 3e3)"]
#write the new template to disc
gtFile <- tempfile()
write.csv(tbl, file = gtFile)
##reload the new template
gt <- gatingTemplate(gtFile, autostart = 1L)
#run the gating
gating(gt, gs)
#hide the gates that are not of interest
toggle.helperGates(gt, gs)
#visualize the gates
library(ggcyto)
autoplot(gs[[1]])

## -----------------------------------------------------------------------------
outFile <- tempfile(fileext = ".xml")
gatingset_to_cytobank(gs, outFile)

## ----eval=FALSE---------------------------------------------------------------
#  outFile <- tempfile(fileext = ".wsp")
#  gatingset_to_flowjo(gs, outFile)

