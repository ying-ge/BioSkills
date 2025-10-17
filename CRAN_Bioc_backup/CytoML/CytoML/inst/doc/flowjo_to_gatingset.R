## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set()

## -----------------------------------------------------------------------------
library(CytoML)
dataDir <- system.file("extdata",package="flowWorkspaceData")
wsfile <- list.files(dataDir, pattern="manual.xml",full=TRUE)
ws <- open_flowjo_xml(wsfile)
ws

## -----------------------------------------------------------------------------
tail(fj_ws_get_sample_groups(ws))

## -----------------------------------------------------------------------------
fj_ws_get_samples(ws, group_id = 5)

## -----------------------------------------------------------------------------
fj_ws_get_keywords(ws, 28)[1:5]

## ----eval=FALSE---------------------------------------------------------------
#  gs <- flowjo_to_gatingset(ws, name = "T-cell")

## ----eval=FALSE---------------------------------------------------------------
#  gs <- flowjo_to_gatingset(ws, name = 4)

## ----eval=FALSE---------------------------------------------------------------
#  gs <- flowjo_to_gatingset(ws, name = 4, path = dataDir)

## -----------------------------------------------------------------------------
gs <- flowjo_to_gatingset(ws, name = 4, execute = FALSE)
gs

## -----------------------------------------------------------------------------
suppressMessages(library(flowWorkspace))
plot(gs)

## -----------------------------------------------------------------------------
gs_pop_get_gate(gs, "CD3+")

## -----------------------------------------------------------------------------
gs_get_compensations(gs)[1]

## -----------------------------------------------------------------------------
gh_get_transformations(gs[[1]], channel = "B710-A")

## -----------------------------------------------------------------------------
head(gs_pop_get_stats(gs, xml = TRUE))

## ----error=TRUE---------------------------------------------------------------
gs_pop_get_data(gs)

## -----------------------------------------------------------------------------
sampleNames(gs)

## -----------------------------------------------------------------------------
gs <- flowjo_to_gatingset(ws, name = 4, execute = FALSE, additional.keys = NULL)
sampleNames(gs)

## -----------------------------------------------------------------------------
gs <- flowjo_to_gatingset(ws, name = 4, execute = FALSE, additional.keys = c("$TOT", "EXPERIMENT NAME"))
sampleNames(gs)


## -----------------------------------------------------------------------------
gs <- flowjo_to_gatingset(ws, name = 4, execute = FALSE, additional.sampleID = TRUE)
sampleNames(gs)

## -----------------------------------------------------------------------------
pData(gs)

## -----------------------------------------------------------------------------
gs <- flowjo_to_gatingset(ws, name = 4, execute = FALSE, additional.keys = NULL, keywords = c("EXPERIMENT NAME", "TUBE NAME"))
pData(gs)

## -----------------------------------------------------------------------------
gs <- flowjo_to_gatingset(ws, name = 4, execute = FALSE, additional.keys = NULL, subset = 1:2)
sampleNames(gs)

## ----eval=FALSE---------------------------------------------------------------
#  gs <- flowjo_to_gatingset(ws, name = 4, execute = FALSE, additional.keys = NULL, subset = c("CytoTrol_CytoTrol_3.fcs"))
#  sampleNames(gs)

## -----------------------------------------------------------------------------
gs <- flowjo_to_gatingset(ws, name = 4, execute = FALSE, additional.keys = NULL
                                                        , subset = `EXPERIMENT NAME` == "C2_Tcell"
                                                        , keywords = c("EXPERIMENT NAME")
                                                        )

pData(gs)

## ----eval=FALSE---------------------------------------------------------------
#  gs <- flowjo_to_gatingset(ws, name = 4, mc.cores = 4)

## ----echo=FALSE---------------------------------------------------------------
is_local <- dir.exists("~/rglab/workspace/CytoML/wsTestSuite")

## ----echo=FALSE, eval=is_local------------------------------------------------
#  path <- "~/rglab/workspace/CytoML/wsTestSuite"
#  thisPath <- file.path(path, "searchRefNode")
#  wsFile <- file.path(thisPath, "2583-Y-MAL067-FJ.xml")
#  ws <- open_flowjo_xml(wsFile)

## ----eval=is_local------------------------------------------------------------
#  gs <- flowjo_to_gatingset(ws, name="Samples", subset = "1379326.fcs", execute = FALSE)
#  nodes <- gs_get_pop_paths(gs)
#  length(nodes)
#  plot(gs, "3+")

## ----eval=is_local------------------------------------------------------------
#  tail(nodes, 10)

## ----eval=is_local------------------------------------------------------------
#  gs <- flowjo_to_gatingset(ws, name="Samples", subset = "1379326.fcs", leaf.bool = F)
#  gs_pop_get_stats(gs)

## ----eval=is_local------------------------------------------------------------
#  recompute(gs)
#  gs_pop_get_stats(gs)

## ----echo=FALSE, eval=is_local------------------------------------------------
#  wsFile <- file.path(path, "bypassfaultynode.xml")
#  ws <- open_flowjo_xml(wsFile)
#  dataDir <- file.path(path,"Cytotrol/NHLBI/Tcell/")

## ----error=TRUE, eval=is_local------------------------------------------------
#  gs <- flowjo_to_gatingset(ws, name = 4, path = dataDir, subset = 1)

## ----eval=is_local------------------------------------------------------------
#  gs <- flowjo_to_gatingset(ws, name = 4, path = dataDir, execute = FALSE)
#  plot(gs)
#  gs_pop_get_gate(gs, "CD4")[[1]]

## ----eval=is_local------------------------------------------------------------
#  gs <- flowjo_to_gatingset(ws, name = 4, path = dataDir, subset = 1, skip_faulty_gate = TRUE)
#  head(gs_pop_get_stats(gs))

## ----echo=FALSE, eval=is_local------------------------------------------------
#   wsFile <- file.path(path, "no-gate.wsp")
#  ws <- open_flowjo_xml(wsFile, sample_names_from = 'sampleNode')

## ----error=TRUE, eval=is_local------------------------------------------------
#  gs <- flowjo_to_gatingset(ws, name = 1, path = dataDir)

## ----eval=is_local------------------------------------------------------------
#  fj_ws_get_samples(ws)

## ----eval=is_local------------------------------------------------------------
#  gs <- flowjo_to_gatingset(ws, name = 1, path = dataDir, include_empty_tree = TRUE)
#  range(gs_cyto_data(gs)[[1]])

## ----echo=FALSE, eval=is_local------------------------------------------------
#  wsFile <- file.path(path, "logicle.wsp")
#  ws <- open_flowjo_xml(wsFile)
#  dataDir <- system.file("extdata", package = "flowCore")

## ----error=TRUE, eval=is_local------------------------------------------------
#  gs <- flowjo_to_gatingset(ws, name = 1, path = dataDir, additional.keys = NULL)

## ----eval=is_local------------------------------------------------------------
#  gs <- flowjo_to_gatingset(ws, name = 1, path = dataDir, additional.keys = NULL, fcs_file_extension = ".B08")
#  sampleNames(gs)

