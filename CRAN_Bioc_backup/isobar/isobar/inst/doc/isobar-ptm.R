### R code from vignette source 'isobar-ptm.Rnw'

###################################################
### code chunk number 1: init
###################################################
  require(ggplot2)
  dir.create(file.path("graphics"), showWarnings = FALSE)


###################################################
### code chunk number 2: load-isobar
###################################################
  library(isobar) ## load the isobar package


###################################################
### code chunk number 3: isobar-ptm.Rnw:81-87 (eval = FALSE)
###################################################
## # Generate PhosphoRS XML input file based on MGF and identification file
## #  massTolerance: fragment ion mass tolerance (in Da)
## #  activationType: CID, HCD, or ETD
## writePhosphoRSInput("phosphors.in.xml",
##                       "identifications.id.csv","peaklist.mgf",
##                        massTolerance=0.5,activationType="CID")


###################################################
### code chunk number 4: isobar-ptm.Rnw:94-98 (eval = FALSE)
###################################################
## # Read PhosphoRS XML output file
## #   simplify: if TRUE, a data.frame is returned, else a list
## #   besthit.only: if TRUE, only the best localization per spectrum is returned
## readPhosphoRSOutput("phosphors.out.xml",simplify=TRUE,besthit.only=TRUE)


###################################################
### code chunk number 5: isobar-ptm.Rnw:104-107 (eval = FALSE)
###################################################
## getPhosphoRSProbabilities("identifications.id.csv","peaklist.mgf",
##                             massTolerance=0.5,activationType="CID",
##                             phosphors.cmd="java -jar phosphoRS.jar")


###################################################
### code chunk number 6: isobar-ptm.Rnw:125-140 (eval = FALSE)
###################################################
## # filterSpectraDeltaScore calls calc.delta.score 
## #   if no column named delta.score is present in the data frame
## # identifications below a min.delta.score are REMOVED
## ib <- readIBSpectra("identifications.id.csv","peaklist.mgf",
##                       annotate.spectra.f=function(...) 
##                         filterSpectraDeltaScore(...,min.delta.score=10))
## 
## 
## # filterSpectraPhosphoRS calls PhosphoRS to caluclate PhosphoRS probabilities
## # identifications below a min.prob (PhosphoRS peptide isoform probability) 
## # are marked to be NOT QUANTIFIED (use.for.quant=FALSE), but not removed
## ib <- readIBSpectra("identifications.id.csv","peaklist.mgf",
##                       annotate.spectra.f=
##                         function(...) filterSpectraPhosphoRS(...,min.prob=0.9,
##                            phosphors.cmd="java -jar PhosphoRS.jar"))


###################################################
### code chunk number 7: isobar-ptm.Rnw:161-165
###################################################
data(ib_phospho)
data(noise.model.hcd)
head(proteinGroup(ib_phospho)@peptideInfo)
10^estimateRatio(ib_phospho,noise.model.hcd,peptide="SPLSPTETFSWPDVR")


###################################################
### code chunk number 8: isobar-ptm.Rnw:170-176
###################################################
pep.n.modif <- unique(apply(fData(ib_phospho)[,c("peptide","modif")],2,cbind))
print(head(pep.n.modif))
estimateRatio(ib_phospho,noise.model.hcd,channel1="114",channel2="115",
                peptide=head(pep.n.modif),combine=FALSE)[,c("lratio","variance",
                                                            "n.spectra","p.value.rat")]



###################################################
### code chunk number 9: ratiodistr
###################################################
suppressPackageStartupMessages(library(distr))
suppressPackageStartupMessages(library(ggplot2))
peptide.ratios <- peptideRatios(ib_phospho,noise.model=noise.model.hcd,
                                  cmbn=matrix(c("114","116"),ncol=1))

lim <- max(abs(peptide.ratios$lratio),na.rm=TRUE)
peptide.distr.cauchy <- fitCauchy(peptide.ratios$lratio)

pseq <- seq(from=-lim,to=lim,length.out=1000)
ggplot() + 
  geom_histogram(aes(x=lratio,y=..density..),data=peptide.ratios,binwidth=0.05,
                 color="darkgreen",fill="white") +
  geom_line(aes(x=x,y=y),color="black",
            data=data.frame(x=pseq,y=d(peptide.distr.cauchy)(pseq)))


###################################################
### code chunk number 10: isobar-ptm.Rnw:233-243
###################################################
peptides <- pep.n.modif[1:5,]

orig.ratio <- estimateRatio(ib_phospho,noise.model.hcd,channel1="114",channel2="115",
                              peptide=peptides,combine=FALSE)[,c("lratio","variance")]

peptides.c <- cbind(peptides,correct.ratio=c(0,-1,1,2,-2))
corr.ratio <- estimateRatio(ib_phospho,noise.model.hcd,channel1="114",channel2="115",
                              peptide=peptides.c,combine=FALSE)[,c("lratio","variance")]

data.frame(peptides.c,orig.ratio,corr.ratio)


###################################################
### code chunk number 11: isobar-ptm.Rnw:264-266 (eval = FALSE)
###################################################
## ptm.info <- getPtmInfoFromPhosphoSitePlus(proteinGroup(ib_phospho),modif="PHOS")
## ptm.info <- getPtmInfoFromNextprot(proteinGroup(ib_phospho))


###################################################
### code chunk number 12: isobar-ptm.Rnw:268-269
###################################################
head(ptm.info)


