## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, results = "markup", message = FALSE, warning = FALSE)

## -----------------------------------------------------------------------------
library(flowWorkspace)
library(CytoML)
acs <- system.file("extdata/cytobank_experiment.acs", package = "CytoML")

## -----------------------------------------------------------------------------
ce <- open_cytobank_experiment(acs)
ce

## -----------------------------------------------------------------------------
sampleNames(ce)
ce_get_panels(ce)
ce_get_compensations(ce)
ce_get_samples(ce)
ce_get_channels(ce)
ce_get_markers(ce)
pData(ce)

## -----------------------------------------------------------------------------
gs <- cytobank_to_gatingset(ce)

## ----eval=FALSE---------------------------------------------------------------
#  xmlfile <- ce$gatingML
#  fcsFiles <- list.files(ce$fcsdir, full.names = TRUE)
#  gs <- cytobank_to_gatingset(xmlfile, fcsFiles)

## -----------------------------------------------------------------------------
library(ggcyto)
## Plot the gates
autoplot(gs[[1]])
# Extract the population statistics
gs_pop_get_count_fast(gs, statType = "count")

