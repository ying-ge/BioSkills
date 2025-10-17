## ----echo=FALSE, results="hide", message=FALSE--------------------------------
require(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)

## ----setup, echo=FALSE, message=FALSE-----------------------------------------
library(batchelor)

## -----------------------------------------------------------------------------
B1 <- matrix(rnorm(10000), ncol=50) # Batch 1 
B2 <- matrix(rnorm(10000), ncol=50) # Batch 2

# Switching easily between batch correction methods.
m.out <- batchCorrect(B1, B2, PARAM=ClassicMnnParam())
f.out <- batchCorrect(B1, B2, PARAM=FastMnnParam(d=20))
r.out <- batchCorrect(B1, B2, PARAM=RescaleParam(pseudo.count=0))

## -----------------------------------------------------------------------------
noCorrect <- function(...) 
# Takes a set of batches and returns them without modification. 
{
   do.call(cbind, list(...)) 
}

## -----------------------------------------------------------------------------
NothingParam <- setClass("NothingParam", contains="BatchelorParam")

## -----------------------------------------------------------------------------
nothing <- NothingParam()
nothing
nothing$some_value <- 1
nothing

## -----------------------------------------------------------------------------
batchCorrect

## -----------------------------------------------------------------------------
setMethod("batchCorrect", "NothingParam", function(..., batch = NULL, 
    restrict=NULL, subset.row = NULL, correct.all = FALSE, 
    assay.type = "logcounts", PARAM) 
{
    batches <- list(...)
    checkBatchConsistency(batches)

    # Pulling out information from the SCE objects.        
    is.sce <- checkIfSCE(batches)
    if (any(is.sce)) {
        batches[is.sce] <- lapply(batches[is.sce], assay, i=assay.type)
    }

    # Subsetting by 'batch', if only one object is supplied. 
    do.split <- length(batches)==1L
    if (do.split) {
        divided <- divideIntoBatches(batches[[1]], batch=batch, restrict=restrict)
        batches <- divided$batches
        restrict <- divided$restricted
    } 

    # Subsetting by row.
    # This is a per-gene "method", so correct.all=TRUE will ignore subset.row.
    # More complex methods will need to handle this differently.
    if (correct.all) {
        subset.row <- NULL
    } else if (!is.null(subset.row)) {
        subset.row <- normalizeSingleBracketSubscript(originals[[1]], subset.row)
        batches <- lapply(batches, "[", i=subset.row, , drop=FALSE)
    }

    # Don't really need to consider restrict!=NULL here, as this function
    # doesn't do anything with the cells anyway.
    output <- do.call(noCorrect, batches)

    # Reordering the output for correctness if it was previously split.
    if (do.split) {
        d.reo <- divided$reorder
        output <- output[,d.reo,drop=FALSE]
    }

    ncells.per.batch <- vapply(batches, FUN=ncol, FUN.VALUE=0L)
    batch.names <- names(batches)
    if (is.null(batch.names)) {
        batch.names <- seq_along(batches)
    }
    
    SingleCellExperiment(list(corrected=output), 
        colData=DataFrame(batch=rep(batch.names, ncells.per.batch)))
})

## -----------------------------------------------------------------------------
n.out <- batchCorrect(B1, B2, PARAM=NothingParam())
n.out

## -----------------------------------------------------------------------------
sessionInfo()

