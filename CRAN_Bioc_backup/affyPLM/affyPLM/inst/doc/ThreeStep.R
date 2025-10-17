### R code from vignette source 'ThreeStep.Rnw'

###################################################
### code chunk number 1: ThreeStep.Rnw:36-37
###################################################
library(affyPLM)


###################################################
### code chunk number 2: ThreeStep.Rnw:44-49 (eval = FALSE)
###################################################
## require(affydata)
## data(Dilution)
## ##FIXME: remove the next line
## Dilution = updateObject(Dilution)
## eset <- threestep(Dilution)


###################################################
### code chunk number 3: ThreeStep.Rnw:54-56 (eval = FALSE)
###################################################
## eset <- threestep(Dilution, background.method = "MASIM",
## 	normalize.method="quantile",summary.method="tukey.biweight")


###################################################
### code chunk number 4: ThreeStep.Rnw:62-64 (eval = FALSE)
###################################################
## eset <- threestep(Dilution, background.method = "IdealMM",
## 	normalize="quantile",summary.method="log.2nd.largest")


