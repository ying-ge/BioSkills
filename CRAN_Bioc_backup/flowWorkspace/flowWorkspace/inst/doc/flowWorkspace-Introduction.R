## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, results = "markup", message = FALSE)

## ----parsews,message=FALSE, warning=FALSE-------------------------------------
library(flowWorkspace)
dataDir <- system.file("extdata",package="flowWorkspaceData")
gs_archive <- list.files(dataDir, pattern = "gs_bcell_auto",full = TRUE)
gs <- load_gs(gs_archive)
gs


## ----sampleNames--------------------------------------------------------------
sampleNames(gs)

## ----subset-------------------------------------------------------------------
gs[1]

## ----plotTree-----------------------------------------------------------------
plot(gs, bool = TRUE)

## ----gs_get_pop_paths-path-1--------------------------------------------------
gs_get_pop_paths(gs, path = 2)

## ----gs_get_pop_paths-path-full-----------------------------------------------
gs_get_pop_paths(gs, path = "full")

## ----gs_get_pop_paths-path-auto-----------------------------------------------
nodelist <- gs_get_pop_paths(gs, path = "auto")
nodelist

## ----gh_pop_get_gate----------------------------------------------------------
node <- nodelist[3]
g <- gs_pop_get_gate(gs, node)
g

## ----getStats-----------------------------------------------------------------
gs_pop_get_stats(gs)[1:10,]

## ----autoplot-nodeName--------------------------------------------------------
library(ggcyto)
autoplot(gs, node)

## ----annotate-----------------------------------------------------------------
d <- data.frame(sample=factor(c("sample 1", "sample 2")),treatment=factor(c("sample","control")) )
pd <- pData(gs)
pd <- cbind(pd,d)
pData(gs) <- pd
pData(gs)

## -----------------------------------------------------------------------------
subset(gs, treatment == "control")

## -----------------------------------------------------------------------------
cs <- gs_pop_get_data(gs)
class(cs)
nrow(cs[[1]])

## ----getData-gh---------------------------------------------------------------
cs <- gs_pop_get_data(gs, node)
nrow(cs[[1]])

## ----gh-----------------------------------------------------------------------
gh <- gs[[1]]
gh

## -----------------------------------------------------------------------------
autoplot(gh)

## ----getInd-------------------------------------------------------------------
table(gh_pop_get_indices(gh,node))

## ----create gs----------------------------------------------------------------
library(flowCore)
data(GvHD)
#select raw flow data
fs <- GvHD[1:2]

## ----GatingSet constructor----------------------------------------------------
gs <- GatingSet(fs)

## ----compensate---------------------------------------------------------------
cfile <- system.file("extdata","compdata","compmatrix", package="flowCore")
comp.mat <- read.table(cfile, header=TRUE, skip=2, check.names = FALSE)
## create a compensation object 
comp <- compensation(comp.mat)
#compensate GatingSet
gs <- compensate(gs, comp)

## ----eval=FALSE---------------------------------------------------------------
#  gs <- compensate(gs, comp.list)

## ----user-transformation------------------------------------------------------
require(scales)
trans.func <- asinh
inv.func <- sinh
trans.obj <- trans_new("myAsinh", trans.func, inv.func)

## ----transform-build-in-------------------------------------------------------
trans.obj <- asinhtGml2_trans()
trans.obj

## ----transformerList----------------------------------------------------------
chnls <- colnames(fs)[3:6] 
transList <- transformerList(chnls, trans.obj)

## ----estimateLogicle----------------------------------------------------------
estimateLogicle(gs[[1]], chnls)

## ----transform-gs-------------------------------------------------------------
gs <- transform(gs, transList)
gs_get_pop_paths(gs) 

## ----add-rectGate-------------------------------------------------------------
rg <- rectangleGate("FSC-H"=c(200,400), "SSC-H"=c(250, 400), filterId="rectangle")
nodeID <- gs_pop_add(gs, rg)
nodeID
gs_get_pop_paths(gs)  

## ----add-quadGate-------------------------------------------------------------
qg <- quadGate("FL1-H"= 0.2, "FL2-H"= 0.4)
nodeIDs <- gs_pop_add(gs,qg,parent="rectangle")
nodeIDs 
gs_get_pop_paths(gs)

## ----add-boolGate-------------------------------------------------------------
bg <- booleanFilter(`CD15 FITC-CD45 PE+|CD15 FITC+CD45 PE-`)
bg
nodeID2 <- gs_pop_add(gs,bg,parent="rectangle")
nodeID2
gs_get_pop_paths(gs)

## ----plot-gh------------------------------------------------------------------
plot(gs, bool=TRUE)

## ----recompute----------------------------------------------------------------
recompute(gs)

## ----autoplot-rect------------------------------------------------------------
autoplot(gs,"rectangle") #plot one Gate

## ----autoplot-multiple--------------------------------------------------------
autoplot(gs, gs_pop_get_children(gs[[1]], "rectangle")[1:4])

## ----autoplot-gh-bool---------------------------------------------------------
autoplot(gs[[1]])

## ----getCMAT------------------------------------------------------------------
gh <- gs[[1]]
gh_get_compensations(gh);

## ----getTrans,results='markup'------------------------------------------------
trans <- gh_get_transformations(gh)
names(trans)
trans[[1]]

## ----rm-----------------------------------------------------------------------
Rm('rectangle', gs)
gs_get_pop_paths(gs)

## ----archive,eval=FALSE-------------------------------------------------------
#  tmp <- tempdir()
#  save_gs(gs,path = file.path(tmp,"my_gs"))
#  gs <- load_gs(file.path(tmp,"my_gs"))

## ----clone,eval=FALSE---------------------------------------------------------
#  gs1 <- gs_clone(gs)

## ----copy-tree,eval=FALSE-----------------------------------------------------
#  gs2 <- gs_copy_tree_only(gs)

## ----load_cf------------------------------------------------------------------
files <- list.files(dataDir, "Cyto", full.names = TRUE)
cf <- load_cytoframe_from_fcs(files[1], num_threads = 4)
cf

## ----load_cf_header-----------------------------------------------------------
cfh <- load_cytoframe_from_fcs(files[1], text.only = TRUE)
cfh

## ----dim_cf-------------------------------------------------------------------
dim(cf)

## ----colnames_cf--------------------------------------------------------------
colnames(cf)

## ----exprs_cf-----------------------------------------------------------------
head(exprs(cf))

## ----spill_cf-----------------------------------------------------------------
spillover(cf)

## ----keys_cf------------------------------------------------------------------
head(keyword(cf))

## ----ref1_cf------------------------------------------------------------------
cf1 <- cf # cf is a reference
colnames(cf1)[1]

## ----ref2_cf------------------------------------------------------------------
colnames(cf1)[1] <- "t"
colnames(cf)[1] # The change affects the original cf object

## ----view_cf------------------------------------------------------------------
cf1 <- cf[1:10, 2:3]
dim(cf1)
exprs(cf)[2,3]
exprs(cf1)[2,2] <- 0 # data change affects the orignal cf
exprs(cf)[2,3]

## ----shallow_cf---------------------------------------------------------------
cf1 <- cf[]

## ----deep_cf------------------------------------------------------------------
cf <- load_cytoframe_from_fcs(files[1], num_threads = 4) # starting fresh
cf1 <- realize_view(cf[1:10, 2:3])
dim(cf1)
exprs(cf)[2,3]
exprs(cf1)[2,2] <- 0 # data change no longer affects the original cf
exprs(cf)[2,3]
exprs(cf1)[2,2] # but does affect the separate data of cf1

## ----coerce_cf----------------------------------------------------------------
fr <- cytoframe_to_flowFrame(cf)
class(fr)
cf_back <- flowFrame_to_cytoframe(fr)
class(cf_back)

## ----pnt_cmp------------------------------------------------------------------
identical(cf@pointer, cf_back@pointer) # These point to distinct copies of the data

## ----h5-----------------------------------------------------------------------
tmpfile <- tempfile(fileext = ".h5")
cf_write_h5(cf, tmpfile)
loaded <- load_cytoframe(tmpfile)

## ----load_cs------------------------------------------------------------------
files <- list.files(dataDir, "Cyto", full.names = TRUE)
cs <- load_cytoset_from_fcs(files, num_threads = 4)
cs

## -----------------------------------------------------------------------------
tmp <- tempfile()
save_cytoset(cs, tmp)
cs <- load_cytoset(tmp, backend_readonly = FALSE)

## ----colnames_cs--------------------------------------------------------------
colnames(cs)

## ----subset_cs----------------------------------------------------------------
sub_cs <- cs[1]

## ----extract1_cf--------------------------------------------------------------
sub_fr <- cs[[1]]
exprs(cs[[1]])[2,2]
exprs(sub_fr)[2,2] <- 0 # This WILL affect the original data
exprs(cs[[1]])[2,2]

## ----extract2_cf--------------------------------------------------------------
sub_cf <- cs[[1, returnType = "flowFrame"]]
exprs(cs[[1]])[2,2]
exprs(sub_cf)[2,2] <- 100 # This WILL NOT affect the original data
exprs(cs[[1]])[2,2]

## ----extract3_cf--------------------------------------------------------------
sub_cf <- get_cytoframe_from_cs(cs,1)

