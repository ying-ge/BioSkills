## ----installPackage, eval=FALSE-----------------------------------------------
#  if (!requireNamespace("BiocManager", quietly=TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("CellNOptR")

## ----installPackage2, eval=FALSE----------------------------------------------
#  if (!requireNamespace("remotes", quietly=TRUE))
#  	install.packages("remotes")
#  remotes::install_github("saezlab/CellNOptR")

## ----Ropts, include=FALSE-------------------------------------------
options(width=70)

## ----loadLib, eval=TRUE, message=FALSE, warning=FALSE---------------
library(CellNOptR)

## ----newDir, eval=FALSE---------------------------------------------
#  dir.create("CNOR_analysis")

## ----quickstart, eval=FALSE-----------------------------------------
#  # ---------------------- load the library and get a SIF and MIDAS file
#  library(CellNOptR)
#  
#  # ---------------------- examples are provided in CellNOptR
#  data("ToyModel", package="CellNOptR")
#  data("CNOlistToy", package="CellNOptR")
#  pknmodel = ToyModel
#  cnolist = CNOlist(CNOlistToy)
#  
#  # ---------------------- alternatively you can read your own files:
#  # pknmodel = readSIF("ToyModel.sif")
#  # cnolist = CNOlist("ToyDataMMB.csv")
#  
#  # ---------------------- preprocess the network
#  model <- preprocessing(cnolist, pknmodel)
#  
#  # ---------------------- perform the analysis
#  res <- gaBinaryT1(CNOlist = cnolist,model =  model, verbose=FALSE)
#  
#  # ---------------------- plot the results
#  cutAndPlot(cnolist, model, list(res$bString))

## ----directory, eval=TRUE, include=FALSE----------------------------
cpfile <- dir(system.file("ToyModel", package="CellNOptR"), full=TRUE)
file.copy(from=cpfile, to=getwd(), overwrite=TRUE)

## ----getData, eval=TRUE---------------------------------------------
dataToy <- readMIDAS("ToyDataMMB.csv", verbose=FALSE)
CNOlistToy <- makeCNOlist(dataToy, subfield=FALSE, verbose=FALSE)

## ----getData2, eval=TRUE--------------------------------------------
data(CNOlistToy,package="CellNOptR", verbose=FALSE)

## ----cnolistClass, eval=TRUE, include=FALSE-------------------------
CNOlistToy = CNOlist("ToyDataMMB.csv")

## ----cnolistClass2, eval=FALSE--------------------------------------
#  data(CNOlistToy,package="CellNOptR")
#  CNOlistToy = CNOlist(CNOlistToy)

## ----showCNO, eval=TRUE---------------------------------------------
CNOlistToy

## ----plotCNO, fig.width=7, fig.height=7, fig.cap="Figure 1: CNOlist data shown by plotting function (either *plot* or *plotCNOlist*)"----
plot(CNOlistToy)

## ----ploCNOPDF, eval=FALSE, include=TRUE----------------------------
#  plotCNOlistPDF(CNOlist=CNOlistToy,filename="ToyModelGraph.pdf")

## ----getModel, eval=TRUE--------------------------------------------
pknmodel<-readSIF("ToyPKNMMB.sif")

## ----getModel2, eval=TRUE-------------------------------------------
data(ToyModel,package="CellNOptR")

## ----indices, eval=TRUE---------------------------------------------
checkSignals(CNOlistToy,pknmodel)

## ----plotModel, eval=TRUE, fig.cap="Figure 2: Prior knowledge network (original SIF file visualised by *plotModel*) for the Toy Model example."----
plotModel(pknmodel, CNOlistToy)

## ----NONC, eval=TRUE------------------------------------------------
indicesToy<-indexFinder(CNOlistToy,pknmodel,verbose=TRUE)
ToyNCNOindices<-findNONC(pknmodel,indicesToy,verbose=TRUE)
ToyNCNOcut<-cutNONC(pknmodel,ToyNCNOindices)
indicesToyNCNOcut<-indexFinder(CNOlistToy,ToyNCNOcut)

## ----compress, eval=TRUE--------------------------------------------
ToyNCNOcutComp<-compressModel(ToyNCNOcut,indicesToyNCNOcut)
indicesToyNCNOcutComp<-indexFinder(CNOlistToy,ToyNCNOcutComp)

## ----expand, eval=TRUE----------------------------------------------
model <- expandGates(ToyNCNOcutComp, maxInputsPerGate=3)

## ----expand2, eval=TRUE---------------------------------------------
model <- preprocessing(CNOlistToy, pknmodel, expansion=TRUE,
    compression=TRUE, cutNONC=TRUE, verbose=FALSE)

## ----resError, eval=TRUE--------------------------------------------
resECNOlistToy <- residualError(CNOlistToy)

## ----initbs, eval=TRUE----------------------------------------------
initBstring <- rep(1,length(model$reacID))

## ----optim, eval=TRUE-----------------------------------------------
ToyT1opt <- gaBinaryT1(CNOlist=CNOlistToy, model=model,
    initBstring=initBstring, verbose=FALSE)

## ----resSim,fig.width=7, fig.height=7, fig.cap="Figure 3: Results of the cutAndPlot function on the Toy Model example."----
cutAndPlot(model=model, bStrings=list(ToyT1opt$bString),
    CNOlist=CNOlistToy,plotPDF=TRUE)

## ----plotFit,fig.width=7, fig.height=7, eval=TRUE, fig.cap="Figure 4: Results of the cutAndPlot function on the Toy Model example"----
plotFit(optRes=ToyT1opt)

## ----simFitPDF, include=TRUE, eval=FALSE----------------------------
#  cutAndPlot(
#      model=model,
#      bStrings=list(ToyT1opt$bString),
#      CNOlist=CNOlistToy,
#      plotPDF=TRUE)
#  pdf("evolFitToyT1.pdf")
#  plotFit(optRes=ToyT1opt)
#  dev.off()

## ----plotModelToy1, eval=T,fig.cap="Figure 4a: Processed model. "----
plotModel(model, CNOlistToy, bString=ToyT1opt$bString)

## ----plotModelToy2, eval=TRUE, fig=F, include=TRUE, fig.cap="Figure 4b: The edges are on (black or red) or off (grey or pink) according to the best set of parameters found during the optimisation (the best bit string). To obtain the right hand side model, the mapBack function has been used."----

bs = mapBack(model, pknmodel, ToyT1opt$bString)
plotModel(pknmodel, CNOlistToy, bs, compressed=model$speciesCompressed)

## ----eval=FALSE-----------------------------------------------------
#  writeScaffold(
#      modelComprExpanded=model,
#      optimResT1=ToyT1opt,
#      optimResT2=NA,
#      modelOriginal=pknmodel,
#      CNOlist=CNOlistToy)
#  
#  writeNetwork(
#      modelOriginal=pknmodel,
#      modelComprExpanded=model,
#      optimResT1=ToyT1opt,
#      optimResT2=NA,
#      CNOlist=CNOlistToy)
#  
#  namesFilesToy<-list(
#      dataPlot="ToyModelGraph.pdf",
#      evolFitT1="evolFitToyT1.pdf",
#      evolFitT2=NA,
#      simResultsT1="SimResultsT1_1.pdf",
#      simResultsT2=NA,
#      scaffold="Scaffold.sif",
#      scaffoldDot="Scaffold.dot",
#      tscaffold="TimesScaffold.EA",
#      wscaffold="weightsScaffold.EA",
#      PKN="PKN.sif",
#      PKNdot="PKN.dot",
#      wPKN="TimesPKN.EA",
#      nPKN="nodesPKN.NA")
#  
#  writeReport(
#      modelOriginal=pknmodel,
#      modelOpt=model,
#      optimResT1=ToyT1opt,
#      optimResT2=NA,
#      CNOlist=CNOlistToy,
#      directory="testToy",
#      namesFiles=namesFilesToy,
#      namesData=list(CNOlist="Toy",model="ToyModel"))

## ----eval=TRUE, echo=FALSE------------------------------------------
unlink("testToy",recursive=TRUE)

## ----eval=TRUE------------------------------------------------------
dataToy <- readMIDAS("ToyDataMMB.csv")
CNOlistToy <- makeCNOlist(dataToy,subfield=FALSE)
pknmodel <- readSIF("ToyPKNMMB.sif")

## ----wrap1, eval=FALSE, include=TRUE--------------------------------
#  res <- CNORwrap(
#      paramsList=NA,
#      name="Toy",
#      namesData=list(CNOlist="ToyData",model="ToyModel"),
#      data=CNOlistToy,
#      model=pknmodel)

## ----eraseToyDir, eval=TRUE, include=FALSE--------------------------
unlink("Toy",recursive=TRUE)

## ----wrap2, eval=FALSE----------------------------------------------
#  pList<-defaultParameters(CNOlistToy, pknmodel)
#  #pList$data = CNOlistToy
#  #pList$model = ToyModel
#  res <- CNORwrap(
#      paramsList=pList,
#      name="Toy1Step",
#      namesData=list(CNOlist="ToyData",model="ToyModel"))

## ----eraseData, eval=TRUE, echo=FALSE-------------------------------
unlink("ToyDataMMB.csv")
unlink("ToyPKNMMB.sif")
unlink("Toy1Step",recursive=TRUE)
unlink("ToyModelMMB2.sif")

## ----getDREAM, eval=TRUE--------------------------------------------
#Option 1: copy the SIF and MIDAS files (followed by readMIDAS, makeCNOlist and readSIF)
cpfile<-dir(system.file("DREAMModel",package="CellNOptR"),full=TRUE)
file.copy(from=cpfile,to=getwd(),overwrite=TRUE)
#Option 2: load the CNOlist and model objects
data(CNOlistDREAM,package="CellNOptR")
data(DreamModel,package="CellNOptR")

## ----DREAMAnalysis, eval=FALSE--------------------------------------
#  model = preprocessing(CNOlistDREAM, DreamModel, verbose=FALSE)
#  res = gaBinaryT1(CNOlistDREAM, model, verbose=FALSE, maxTime=10)
#  cutAndPlot(CNOlistDREAM, model, bStrings=list(res$bString), plotPDF=TRUE,
#      plotParams=list(maxrow=25, margin=0.1, width=20, height=20))
#  
#  

## ----t2load, eval=TRUE----------------------------------------------
data(CNOlistToy2,package="CellNOptR")
data(ToyModel2,package="CellNOptR")
pknmodel = ToyModel2
cnolist = CNOlist(CNOlistToy2)

## ----t2plotCNOlist, fig.width=7, eval=TRUE--------------------------
plot(cnolist)

## ----t2Opt1, eval=TRUE----------------------------------------------
model = preprocessing(cnolist, pknmodel, verbose=FALSE)
T1opt <- gaBinaryT1(cnolist, model, stallGenMax=10, maxTime=60, verbose=FALSE)

## ----t2OptT1plot, eval=TRUE, fig.width=7----------------------------
sim <- cutAndPlot(model=model, bStrings=list(T1opt$bString),
    CNOlist=cnolist, plotPDF=TRUE)

## ----eval=TRUE, fig.width=6 , fig.height=6--------------------------
plotFit(optRes=T1opt)

## ----t2OptT2, eval=TRUE---------------------------------------------
T2opt<-gaBinaryTN(cnolist, model,  bStrings=list(T1opt$bString),
    stallGenMax=10, maxTime=60, verbose=FALSE)

## ----resSimT2, eval=FALSE-------------------------------------------
#  res <- cutAndPlot(
#      model=model,
#      bStrings=list(T1opt$bString, T2opt$bString),
#      CNOlist=cnolist,
#      plotPDF=TRUE, plotParams=list(cex=0.8, cmap_scale=0.5, margin=0.2))

## ----eval=FALSE-----------------------------------------------------
#  # ---------------------- load the library and get a SIF and MIDAS file
#  library(CellNOptR)
#  library(stringr)
#  #
#  # ---------------------- examples are provided in CellNOptR
#  data("ToyModel", package="CellNOptR")
#  data("CNOlistToy", package="CellNOptR")
#  pknmodel = ToyModel
#  cnolist = CNOlist(CNOlistToy)
#  #
#  # ---------------------- alternatively you can read your own files:
#  # pknmodel = readSIF("ToyModel.sif")
#  # cnolist = CNOlist("ToyDataMMB.csv")
#  #
#  # ---------------------- preprocess the network
#  model = preprocessing(cnolist, pknmodel)
#  #
#  # ---------------------- perform the analysis
#  cplexPath = "path/to/cplex"
#  resILP = ilpBinaryT1(cnolist = cnolist, model = model,
#                       numSolutions = 3, relGap = 0.05,
#                       cplexPath = cplexPath) # asking to retrieve 3
#                                              # equivalent solutions with
#                                              # tolerance gap relGap=0.05
#  #
#  # ---------------------- plot the results (optimized models + fits)
#  cutAndPlot(CNOlist = cnolist, model = model, bStrings = list(resILP$bitstringILP[[1]]))
#  plotModel(model, cnolist, bString=resILP$bitstringILP[[1]])
#  
#  cutAndPlot(CNOlist = cnolist, model = model, bStrings = list(resILP$bitstringILP[[2]]))
#  plotModel(model, cnolist, bString=resILP$bitstringILP[[2]])
#  
#  cutAndPlot(CNOlist = cnolist, model = model, bStrings = list(resILP$bitstringILP[[3]]))
#  plotModel(model, cnolist, bString=resILP$bitstringILP[[3]])

## ----eval=FALSE-----------------------------------------------------
#  # ---------------------- load the library and get a SIF and MIDAS file
#  library(CellNOptR)
#  #
#  # ---------------------- examples are provided in CellNOptR
#  data(PKN_ToyPB, package="CellNOptR")
#  data(CNOlist_ToyPB, package="CellNOptR")
#  #
#  # ---------------------- preprocess the network
#  model = preprocessing(data = cnodata,model = pknmodel,compression = T,expansion = T)
#  plotModel(model,cnodata)
#  #
#  #------------------------ original CNOlist contains many timepoints, we use only a subset
#  plot(cnodata)
#  selectedTime = c(0,10)
#  cnodata_prep = cutCNOlist(cnodata ,model = model,
#                            cutTimeIndices = which(!cnodata@timepoints %in% selectedTime))
#  plot(cnodata_prep)
#  #
#  # ---------------------- perform the analysis
#  opt = gaBinaryT1(CNOlist = cnodata_prep,model = model)
#  #
#  # ---------------------- 10-fold cross-validation procedure using T1 data
#  library(doParallel)
#  doParallel::registerDoParallel(cores=3)
#  system.time({R1=crossvalidateBoolean(CNOlist = cnodata_prep, model = model,
#                                       type="datapoint", nfolds=10, parallel = TRUE)})
#  system.time({R2=crossvalidateBoolean(CNOlist = cnodata_prep, model=model,
#                                       type="experiment", nfolds=10, parallel = TRUE)})
#  system.time({R3=crossvalidateBoolean(CNOlist = cnodata_prep, model=model,
#                                       type="observable", nfolds=10, parallel = TRUE)})
#  system.time({R4=crossvalidateBoolean(CNOlist = cnodata_prep, model=model,
#                                       type="datapoint", nfolds=10, parallel = FALSE)})
#  system.time({R5=crossvalidateBoolean(CNOlist = cnodata_prep, model=model,
#                                       type="experiment", nfolds=10, parallel = FALSE)})
#  system.time({R6=crossvalidateBoolean(CNOlist = cnodata_prep, model=model,
#                                       type="observable", nfolds=10, parallel = FALSE)})

