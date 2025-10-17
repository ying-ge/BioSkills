## ----eval=FALSE---------------------------------------------------------------
#  gs_split_by_tree(x)
#  gs_check_redundant_nodes(x)
#  gs_remove_redundant_nodes(x,toRemove)
#  gs_remove_redundant_channels(gs, ...)

## ----echo=FALSE, message=FALSE, results='hide'--------------------------------
library(flowWorkspace)
flowDataPath <- system.file("extdata", package = "flowWorkspaceData")
gs <- load_gs(file.path(flowDataPath,"gs_manual"))
gs1 <- gs_clone(gs)
sampleNames(gs1) <- "1.fcs"

# simply the tree
nodes <- gs_get_pop_paths(gs1)
for(toRm in nodes[grepl("CCR", nodes)])
  gs_pop_remove(gs1, toRm)

# remove two terminal nodes
gs2 <- gs_clone(gs1)
sampleNames(gs2) <- "2.fcs"
gs_pop_remove(gs2, "DPT")
gs_pop_remove(gs2, "DNT")

# remove singlets gate
gs3 <- gs_clone(gs2)
gs_pop_remove(gs3, "singlets")
gs_pop_add(gs3, gs_pop_get_gate(gs2, "CD3+"), parent = "not debris")
for(tsub in c("CD4", "CD8"))
  {
    gs_pop_add(gs3, gs_pop_get_gate(gs2, tsub), parent = "CD3+")
    for(toAdd in gs_pop_get_children(gs2, tsub))
    {
        thisParent <- gs_pop_get_parent(gs2[[1]], toAdd,path="auto")
        gs_pop_add(gs3, gs_pop_get_gate(gs2, toAdd), parent = thisParent) 
    }
  }
sampleNames(gs3) <- "3.fcs"

# spin the branch to make it isomorphic
gs4 <- gs_clone(gs3)
# rm cd4 branch first
gs_pop_remove(gs4, "CD4")
# add it back
gs_pop_add(gs4, gs_pop_get_gate(gs3, "CD4"), parent = "CD3+")
# add all the chilren back
for(toAdd in gs_pop_get_children(gs3, "CD4"))
{
    thisParent <- gs_pop_get_parent(gs3[[1]], toAdd)
    gs_pop_add(gs4, gs_pop_get_gate(gs3, toAdd), parent = thisParent)
}
sampleNames(gs4) <- "4.fcs"

gs5 <- gs_clone(gs4)
# add another redundant node
gs_pop_add(gs5, gs_pop_get_gate(gs, "CD4/CCR7+ 45RA+")[[1]], parent = "CD4")
gs_pop_add(gs5, gs_pop_get_gate(gs, "CD4/CCR7+ 45RA-")[[1]], parent = "CD4")
sampleNames(gs5) <- "5.fcs"

library(knitr)
opts_chunk$set(fig.show = 'hold', fig.width = 4, fig.height = 4, results= 'asis')


## ----echo=FALSE---------------------------------------------------------------
plot(gs1)
plot(gs2)

## ----echo=FALSE---------------------------------------------------------------
plot(gs2)
plot(gs3)

## -----------------------------------------------------------------------------
invisible(gs_pop_set_visibility(gs2, "singlets", FALSE))
plot(gs2)
plot(gs3)

## ----results='hold'-----------------------------------------------------------
gs_get_pop_paths(gs2)[5]
gs_get_pop_paths(gs3)[5]

## ----results='hold'-----------------------------------------------------------
gs_get_pop_paths(gs2, path = "auto")[5]
gs_get_pop_paths(gs3, path = "auto")[5]

## ----echo=FALSE---------------------------------------------------------------
#restore gs2
invisible(gs_pop_set_visibility(gs2, "singlets", TRUE))

## ----echo=FALSE---------------------------------------------------------------
plot(gs3)
plot(gs4)

## -----------------------------------------------------------------------------
gslist <- list(gs1, gs2, gs3, gs4, gs5)
gs_groups <- gs_split_by_tree(gslist)
length(gs_groups)

## ----error=TRUE---------------------------------------------------------------
res <- try(gs_check_redundant_nodes(gs_groups), silent = TRUE)
print(res[[1]])

## -----------------------------------------------------------------------------
for(gp in gs_groups)
  plot(gp[[1]])

## -----------------------------------------------------------------------------
for(i in c(2,4))
  for(gs in gs_groups[[i]])
    invisible(gs_pop_set_visibility(gs, "singlets", FALSE))

## -----------------------------------------------------------------------------
toRm <- gs_check_redundant_nodes(gs_groups)
toRm

## ----message=FALSE------------------------------------------------------------
gs_plot_diff_tree(gs_groups)

## ----results='hide'-----------------------------------------------------------
gs_remove_redundant_nodes(gs_groups, toRm)

## -----------------------------------------------------------------------------
GatingSetList(gslist)

## -----------------------------------------------------------------------------
gs_remove_redundant_channels(gs1)

