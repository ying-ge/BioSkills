### R code from vignette source 'isobar.Rnw'

###################################################
### code chunk number 1: init
###################################################
  require(distr)
  require(ggplot2)
  dir.create(file.path("graphics"), showWarnings = FALSE)


###################################################
### code chunk number 2: isobar.Rnw:55-58 (eval = FALSE)
###################################################
## if (!requireNamespace("BiocManager", quietly=TRUE))
##     install.packages("BiocManager")
## BiocManager::install("isobar")


###################################################
### code chunk number 3: isobar.Rnw:67-69 (eval = FALSE)
###################################################
## library(devtools)
## install_github("fbreitwieser/isobar")


###################################################
### code chunk number 4: load-isobar
###################################################
  library(isobar) ## load the isobar package


###################################################
### code chunk number 5: isobar.Rnw:79-80 (eval = FALSE)
###################################################
## packageVersion("isobar")


###################################################
### code chunk number 6: isobar.Rnw:83-84 (eval = FALSE)
###################################################
## packageDescription("isobar")$Version


###################################################
### code chunk number 7: isobar.Rnw:116-126 (eval = FALSE)
###################################################
##   ## generating IBSpectra object from ID.CSV and MGF
##   ib <- readIBSpectra("iTRAQ4plexSpectra",list.files(pattern=".id.csv"),
##           list.files(pattern=".mgf"))
## 
##   ## write in tabular IBSPECTRA.CSV format to file
##   write.table(as.data.frame(ib),sep="\t",row.names=F,
##           file="myexperiment.ibspectra.csv")
## 
##   ## generate from saved IBSPECTRA.CSV - MGF does not have to be supplied
##   ib.2 <- readIBSpectra("iTRAQ4plexSpectra","myexperiment.ibspectra.csv")


###################################################
### code chunk number 8: load-ibspiked
###################################################
  data(ibspiked_set1)
  ceru.human <- protein.g(proteinGroup(ibspiked_set1),"CERU_HUMAN")
  ceru.rat <- protein.g(proteinGroup(ibspiked_set1),"CERU_RAT")
  ceru.mouse <- protein.g(proteinGroup(ibspiked_set1),"CERU_MOUSE")
  ceru.proteins <- c(ceru.human,ceru.rat,ceru.mouse)


###################################################
### code chunk number 9: isobar.Rnw:233-238 (eval = FALSE)
###################################################
## readIBSpectra(...,
##               mapping.file="mapping2.csv",
##               mapping=c(identification.spectrum="id-spectrum-ms2",
##                         quantification.spectrum="quant-spectrum-ms3")
##               )


###################################################
### code chunk number 10: isobar.Rnw:251-258 (eval = FALSE)
###################################################
## readIBSpectra("TMT6plexSpectra",
##               id.file="cid.identifications.csv",
##               peaklist.file="hcd.peaklist.mgf",
##               id.file.domap="hcd.identifications.csv",
##               mapping.file="mapping.csv",
##               ...
##               )


###################################################
### code chunk number 11: isobar.Rnw:275-281 (eval = FALSE)
###################################################
## precursors <- read.delim("file1_precursors.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
## idfile <- read.delim("file1.id.csv",header=TRUE,sep="\t",stringsAsFactors=FALSE)
## 
## idfile <- merge(idfile,precursors,
##                 by.x=c("fixed_spectrum","scans.from"),
##                 by.y=c("fixed_spectrum","SCAN"),all.x=TRUE) # FIXME see /work/analysis/gwinter/raw-spectra/analysis.R


###################################################
### code chunk number 12: isobar.Rnw:300-302 (eval = FALSE)
###################################################
##   as(ibspectra,"MSnSet")
##   as(msnset,"IBSpectra")


###################################################
### code chunk number 13: isobar.Rnw:336-338 (eval = FALSE)
###################################################
##   as(ibspectra,"MSnSet")
##   as(msnset,"IBSpectra")


###################################################
### code chunk number 14: show-reporterMasses
###################################################
sprintf("%.4f",reporterTagMasses(ibspiked_set1))  ## expected masses
mass <- assayData(ibspiked_set1)[["mass"]] ## observerd masses
apply(mass,2,function(x) sprintf("%.4f",quantile(x,na.rm=TRUE,probs=c(0.025,0.975))))


###################################################
### code chunk number 15: print-reporterMasses
###################################################
  print(reporterMassPrecision(ibspiked_set1))


###################################################
### code chunk number 16: normalize-ibspiked
###################################################
ib.old <- ibspiked_set1
ibspiked_set1 <- correctIsotopeImpurities(ibspiked_set1)
ibspiked_set1 <- normalize(ibspiked_set1)


###################################################
### code chunk number 17: fig-maplot-normalize1
###################################################
png("graphics/fig_maplot_normalize.png",width=6.6,height=3,units="in",res=300,pointsize=8)


###################################################
### code chunk number 18: fig-maplot-normalize2
###################################################
par(mfrow=c(1,2))
maplot(ib.old,channel1="114",channel2="117",ylim=c(0.5,2),
       main="before normalization")
abline(h=1,col="red",lwd=2)
maplot(ibspiked_set1,channel1="114",channel2="117",ylim=c(0.5,2),
       main="after normalization")
abline(h=1,col="red",lwd=2)


###################################################
### code chunk number 19: fig-maplot-normalize3
###################################################
dev.off()


###################################################
### code chunk number 20: calc-noisemodel
###################################################
ib.background <- subsetIBSpectra(ibspiked_set1,protein=ceru.proteins,direction="exclude")
noise.model <- NoiseModel(ib.background)


###################################################
### code chunk number 21: ceru-noisemodel
###################################################
ib.ceru <- subsetIBSpectra(ibspiked_set1,protein=ceru.proteins,
                            direction="include",
                            specificity="reporter-specific")
nm.ceru <- NoiseModel(ib.ceru,one.to.one=FALSE,pool=TRUE)


###################################################
### code chunk number 22: maplot-noisemodel1
###################################################
png("graphics/fig_maplot_noisemodel.png",width=3.3,height=3,units="in",res=300,pointsize=8)


###################################################
### code chunk number 23: maplot-noisemodel2
###################################################
maplot(ib.background,noise.model=c(noise.model,nm.ceru),
       channel1="114",channel2="115",ylim=c(0.2,5),
       main="95% CI noise model")


###################################################
### code chunk number 24: maplot-noisemodel3
###################################################
dev.off()


###################################################
### code chunk number 25: estimateRatio
###################################################
## Calculate ratio based on all spectra of peptides specific 
##  to CERU_HUMAN, CERU_RAT or CERU_MOUSE. Returns a named
##  numeric vector.
10^estimateRatio(ibspiked_set1,noise.model,
                 channel1="114",channel2="115",
                 protein=ceru.proteins)['lratio']

## If argument 'combine=FALSE', estimateRatio returns a data.frame 
##  with one row per protein
10^estimateRatio(ibspiked_set1,noise.model,
                 channel1="114",channel2="115",
                 protein=ceru.proteins,combine=FALSE)[,'lratio']
## spiked material channel 115 vs 114: 
##                 CERU_HUMAN (P00450): 1:1
##                 CERU_RAT   (P13635): 2:1  = 2
##                 CERU_MOUSE (Q61147): 5:10 = 0.5

## Peptides shared between rat and mouse
pep.shared <- peptides(proteinGroup(ibspiked_set1),
                       c(ceru.rat,ceru.mouse),set="intersect",
                       columns=c('peptide','n.shared.groups'))
## remove those which are shared with other proteins
pep.shared <- pep.shared$peptide[pep.shared$n.shared.groups==2]

## calculate ratio: it is between the rat and mouse ratios
10^estimateRatio(ibspiked_set1,noise.model,
                 channel1="114",channel2="115",
                 peptide=pep.shared)['lratio']


###################################################
### code chunk number 26: protein-ratios
###################################################
protein.ratios <- proteinRatios(ibspiked_set1,noise.model,cl=c("1","0","0","0"))

## defined class 114 and 115 as class 'T', 116 and 117 as class 'C'
classLabels(ibspiked_set1) <- c("T","T","C","C")

proteinRatios(ibspiked_set1,noise.model,protein=ceru.proteins,
              cl=classLabels(ibspiked_set1),combn.method="interclass",
              summarize=T)[,c("ac","lratio","variance")]


###################################################
### code chunk number 27: ratiodistribution
###################################################
#protein.ratios <- proteinRatios(ibspiked_set1,noise.model)
protein.ratiodistr.wn <- fitWeightedNorm(protein.ratios[,'lratio'],
                                         weights=1/protein.ratios[,'variance'])
protein.ratiodistr.cauchy <- fitCauchy(protein.ratios[,"lratio"])


###################################################
### code chunk number 28: plot-ratiodistr
###################################################
library(distr) # required library
limits=seq(from=-0.5,to=0.5,by=0.001)
curve.wn <- data.frame(x=limits,y=d(protein.ratiodistr.wn)(limits))
curve.cauchy<-data.frame(x=limits,y=d(protein.ratiodistr.cauchy)(limits))

g <- ggplot(data.frame(protein.ratios),aes(x=lratio)) +
  geom_histogram(colour = "darkgreen", fill = "white",aes(y=..density..),
                 binwidth=0.02) + geom_rug() +
  geom_line(data=curve.wn,aes(x=x,y=y),colour="blue") +
  geom_line(data=curve.cauchy,aes(x=x,y=y),colour="red")
print(g)


###################################################
### code chunk number 29: significant-ratios
###################################################
rat.list <-
  estimateRatio(ibspiked_set1,noise.model=noise.model,channel1="114",channel2="115",
                protein=reporterProteins(proteinGroup(ibspiked_set1)),combine=F,
                ratiodistr=protein.ratiodistr.cauchy)
rat.list[rat.list[,"is.significant"]==1,]


###################################################
### code chunk number 30: shared-pep
###################################################
## peptides shared between CERU_RAT and CERU_MOUSE have been computed before
pep.shared
## peptides specific to CERU_RAT
pep.rat <- peptides(proteinGroup(ibspiked_set1),protein=ceru.rat,
                    specificity="reporter-specific")

## create an IBSpectra object with only CERU_RAT and shared peptides
ib.subset <- subsetIBSpectra(ibspiked_set1,
                             peptide=c(pep.rat,pep.shared),direction="include")

## calculate shared ratios
sr <- shared.ratios(ib.subset,noise.model,
                    channel1="114",channel2="117",
                    ratiodistr=protein.ratiodistr.cauchy)
sr



###################################################
### code chunk number 31: isobar.Rnw:631-634
###################################################
  ## plot significantly different protein groups where 90% CI does not overlap
  ## CERU_MOUSE and CERU_RAT is detected, as expected.
  shared.ratios.sign(sr,z.shared=1.282)


###################################################
### code chunk number 32: isobar.Rnw:658-660 (eval = FALSE)
###################################################
## create.reports(type="iTRAQ4plexSpectra",
##                identifications="my.id.csv",peaklist="my.mgf")


###################################################
### code chunk number 33: isobar.Rnw:684-687 (eval = FALSE)
###################################################
## ## execute to find the path and file location in your installation.
## system.file("report",package="isobar") ## path
## list.files(system.file("report",package="isobar")) ## files


###################################################
### code chunk number 34: isobar.Rnw:761-772
###################################################
#p <- readLines(system.file("report","properties.R",package="isobar"))
#p <- sub("(#.*)","\\\\textcolor{comment}{\\1}",p)
#for (i in c("TRUE","FALSE","NULL"))
#  p <- sub(i,paste("\\\\textcolor{special}{",i,"}",sep=""),p)
#p <- sub('(".*?")',"\\\\textcolor{string}{\\1}",p)
#p <- gsub("#","\\#",p,fixed=TRUE)
#p <- gsub("_","\\_",p,fixed=TRUE)
#cat("\\begin{lstlisting}")
#cat(p,sep="\n")
#cat("\\end{lstlisting}")
cat(sprintf("\\lstinputlisting{%s}\n",system.file("report","properties.R",package="isobar")))


###################################################
### code chunk number 35: isobar.Rnw:793-796 (eval = FALSE)
###################################################
## ## execute to find the path and file location in your installation.
## system.file("pl",package="isobar") ## path
## list.files(system.file("pl",package="isobar")) ## files


###################################################
### code chunk number 36: isobar.Rnw:804-810 (eval = FALSE)
###################################################
## ## execute on your system
## system(paste("perl",system.file("pl","mascotParser2.pl",package="isobar"),
##             "--help"))
## 
## print(paste("perl",system.file("pl","pidresParser2.pl",package="isobar"),
##             "--help"))


###################################################
### code chunk number 37: isobar.Rnw:818-819
###################################################
  toLatex(sessionInfo())


