### R code from vignette source 'AffyExtensions.Rnw'

###################################################
### code chunk number 1: loadData
###################################################
library(affyPLM)
require(affydata)
data(Dilution)   # an example dataset provided by the affydata package
##FIXME: drop this after Dilution is updated
Dilution = updateObject(Dilution)
options(width=36)


###################################################
### code chunk number 2: defaultModel
###################################################
Pset <- fitPLM(Dilution)


###################################################
### code chunk number 3: accessors
###################################################
coefs(Pset)[1:5,]
se(Pset)[1:5,]


###################################################
### code chunk number 4: seeDefaultOutput
###################################################
verify.output.param()


###################################################
### code chunk number 5: noResidNoWeights
###################################################
Pset <- fitPLM(Dilution,output.param=list(residuals=FALSE,weights=FALSE))


###################################################
### code chunk number 6: noRobustness
###################################################
Pset <- fitPLM(Dilution,model.param=list(max.its=0))


###################################################
### code chunk number 7: treatmenteffect
###################################################
Pset <- fitPLM(Dilution,  MM ~ -1 + liver + scanner + probes,subset = geneNames(Dilution)[1:100])


###################################################
### code chunk number 8: treatmenteffectexamine
###################################################
coefs(Pset)[1,]


###################################################
### code chunk number 9: treatmenteffectexamine2
###################################################
coefs.probe(Pset)[1]


###################################################
### code chunk number 10: treatmenteffectcovariate
###################################################
logliver <- log2(c(20,20,10,10))
Pset <- fitPLM(Dilution,model=PM~-1+probes+logliver+scanner, variable.type=c(logliver="covariate"),subset = geneNames(Dilution)[1:100])
coefs(Pset)[1,]


###################################################
### code chunk number 11: MMcovariate
###################################################
Pset <- fitPLM(Dilution,  PM ~ MM + samples + probes,subset = geneNames(Dilution)[1:100])


###################################################
### code chunk number 12: MMcovariateexamine
###################################################
coefs(Pset)[1,]
coefs.const(Pset)[1,]
coefs.probe(Pset)[1]


###################################################
### code chunk number 13: probetype
###################################################
Pset <- fitPLM(Dilution,  PMMM ~ liver + probe.type + probes,subset = geneNames(Dilution)[1:100])


###################################################
### code chunk number 14: probetypeexamine
###################################################
coefs(Pset)[1,]
coefs.const(Pset)[1,]
coefs.probe(Pset)[1]


###################################################
### code chunk number 15: probesInTreatment
###################################################
Pset <- fitPLM(Dilution,  PM ~ -1 + liver + liver:probes,subset = geneNames(Dilution)[1:100])


###################################################
### code chunk number 16: probesInTreatmentexamine
###################################################
coefs.probe(Pset)[1]


###################################################
### code chunk number 17: probesInProbetype
###################################################
Pset <- fitPLM(Dilution,   PMMM ~ -1 + liver + probe.type:probes,subset = geneNames(Dilution)[1:100])
coefs.probe(Pset)[1]


###################################################
### code chunk number 18: probesInProbetypeInTreatment
###################################################
Pset <- fitPLM(Dilution,   PMMM ~ -1 + liver + liver:probe.type:probes,subset = geneNames(Dilution)[1:100])
coefs.probe(Pset)[1]


###################################################
### code chunk number 19: constraintExample
###################################################
data(Dilution)
##FIXME: remove next line
Dilution = updateObject(Dilution)
Pset <- fitPLM(Dilution, model = PM ~ probes + samples,constraint.type=c(samples="contr.sum"),subset = geneNames(Dilution)[1:100])
coefs.const(Pset)[1:2]
coefs(Pset)[1:2,]


