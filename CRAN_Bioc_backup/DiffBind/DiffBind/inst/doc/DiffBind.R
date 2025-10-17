### R code from vignette source 'DiffBind.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: DiffBind.Rnw:225-229
###################################################
tmp <-  tempfile(as.character(Sys.getpid()))
pdf(tmp)
savewarn <- options("warn")
options(warn=-1)


###################################################
### code chunk number 3: DiffBind.Rnw:233-234
###################################################
library(DiffBind)


###################################################
### code chunk number 4: DiffBind.Rnw:236-237 (eval = FALSE)
###################################################
## setwd(system.file('extra',package='DiffBind'))


###################################################
### code chunk number 5: DiffBind.Rnw:245-252 (eval = FALSE)
###################################################
## tmpdir <- tempdir()
## url <- 'https://content.cruk.cam.ac.uk/bioinformatics/software/DiffBind/DiffBind_vignette_data.tar.gz'
## file <- basename(url)
## options(timeout=600)
## download.file(url, file.path(tmpdir,file))
## untar(file.path(tmpdir,file), exdir = tmpdir )
## setwd(file.path(tmpdir,"DiffBind_Vignette"))


###################################################
### code chunk number 6: DiffBind.Rnw:258-259 (eval = FALSE)
###################################################
## tamoxifen <- dba.analyze("tamoxifen.csv")


###################################################
### code chunk number 7: DiffBind.Rnw:265-266 (eval = FALSE)
###################################################
## tamoxifen.DB <- dba.report(tamoxifen)


###################################################
### code chunk number 8: DiffBind.Rnw:277-283 (eval = FALSE)
###################################################
## tamoxifen <- dba(sampleSheet="tamoxifen.csv") %>%
##   dba.blacklist() %>%
##   dba.count()     %>%
##   dba.normalize() %>%
##   dba.contrast()  %>%
##   dba.analyze()


###################################################
### code chunk number 9: sampSheet
###################################################
samples <- read.csv(file.path(system.file("extra", package="DiffBind"),
                              "tamoxifen.csv"))
names(samples)
samples


###################################################
### code chunk number 10: dbaConstruct
###################################################
tamoxifen <- dba(sampleSheet="tamoxifen.csv",
                 dir=system.file("extra", package="DiffBind"))


###################################################
### code chunk number 11: dbaConstructDF (eval = FALSE)
###################################################
## tamoxifen <- dba(sampleSheet=samples)


###################################################
### code chunk number 12: DiffBind.Rnw:324-325
###################################################
tamoxifen


###################################################
### code chunk number 13: tamox_occ_corhm
###################################################
plot(tamoxifen)


###################################################
### code chunk number 14: DiffBind.Rnw:377-378 (eval = FALSE)
###################################################
## tamoxifen <- dba.count(tamoxifen) 


###################################################
### code chunk number 15: DiffBind.Rnw:380-381
###################################################
data(tamoxifen_counts)


###################################################
### code chunk number 16: DiffBind.Rnw:389-390
###################################################
tamoxifen


###################################################
### code chunk number 17: eff_lib_size
###################################################
info <- dba.show(tamoxifen)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP, 
                  PeakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID
libsizes


###################################################
### code chunk number 18: tamox_aff_corhm
###################################################
plot(tamoxifen)


###################################################
### code chunk number 19: normalize
###################################################
tamoxifen <- dba.normalize(tamoxifen)


###################################################
### code chunk number 20: show_norm
###################################################
norm <- dba.normalize(tamoxifen, bRetrieve=TRUE)
norm


###################################################
### code chunk number 21: norm_facs
###################################################
normlibs <- cbind(FullLibSize=norm$lib.sizes, NormFacs=norm$norm.factors,
                  NormLibSize=round(norm$lib.sizes/norm$norm.factors))
rownames(normlibs) <- info$ID
normlibs


###################################################
### code chunk number 22: DiffBind.Rnw:500-503
###################################################
tamoxifen <- dba.contrast(tamoxifen, 
                          reorderMeta=list(Condition="Responsive"))
tamoxifen


###################################################
### code chunk number 23: DiffBind.Rnw:535-537
###################################################
tamoxifen <- dba.analyze(tamoxifen)
dba.show(tamoxifen, bContrasts=TRUE)


###################################################
### code chunk number 24: tamox_sdb_corhm
###################################################
plot(tamoxifen, contrast=1)


###################################################
### code chunk number 25: DiffBind.Rnw:587-588
###################################################
tamoxifen.DB <- dba.report(tamoxifen)


###################################################
### code chunk number 26: DiffBind.Rnw:594-595
###################################################
tamoxifen.DB


###################################################
### code chunk number 27: DiffBind.Rnw:615-617
###################################################
sum(tamoxifen.DB$Fold>0)
sum(tamoxifen.DB$Fold<0)


###################################################
### code chunk number 28: tamox_sdb_venn
###################################################
dba.plotVenn(tamoxifen, contrast=1, bDB=TRUE,
             bGain=TRUE, bLoss=TRUE, bAll=FALSE)


###################################################
### code chunk number 29: tamox_aff_pca
###################################################
dba.plotPCA(tamoxifen,DBA_TISSUE,label=DBA_CONDITION)


###################################################
### code chunk number 30: tamox_sdb_pca
###################################################
dba.plotPCA(tamoxifen, contrast=1, label=DBA_TISSUE)


###################################################
### code chunk number 31: tamox_sdb_ma
###################################################
dba.plotMA(tamoxifen)


###################################################
### code chunk number 32: tamox_sdb_volcano
###################################################
dba.plotVolcano(tamoxifen)


###################################################
### code chunk number 33: DiffBind.Rnw:786-788
###################################################
sum(tamoxifen.DB$Fold<0)
sum(tamoxifen.DB$Fold>0)


###################################################
### code chunk number 34: tamox_sdb_box
###################################################
pvals <- dba.plotBox(tamoxifen)


###################################################
### code chunk number 35: DiffBind.Rnw:827-828
###################################################
pvals


###################################################
### code chunk number 36: DiffBind.Rnw:847-848
###################################################
corvals <- dba.plotHeatmap(tamoxifen)


###################################################
### code chunk number 37: tamox_sdb_hm
###################################################
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)
readscores <- dba.plotHeatmap(tamoxifen, contrast=1, correlations=FALSE,
                              scale="row", colScheme = hmap)


###################################################
### code chunk number 38: DiffBind.Rnw:957-959 (eval = FALSE)
###################################################
## profiles <- dba.plotProfile(tamoxifen)
## dba.plotProfile(profiles)


###################################################
### code chunk number 39: DiffBind.Rnw:993-995 (eval = FALSE)
###################################################
## profiles <- dba.plotProfile(tamoxifen,merge=c(DBA_TISSUE, DBA_REPLICATE))
## dba.plotProfile(profiles)


###################################################
### code chunk number 40: DiffBind.Rnw:1018-1028 (eval = FALSE)
###################################################
## mask.MCF7 <- tamoxifen$masks$MCF7
## mask.Resistant  <- tamoxifen$masks$Resistant
## mask.Responsive <- tamoxifen$masks$Responsive
## profiles <- dba.plotProfile(tamoxifen,
##                            samples=list(MCF7_Resistant=
##                                           mask.MCF7 & mask.Resistant,
##                                         MCF7_Responsive=
##                                           mask.MCF7 & mask.Responsive),
##                            merge=NULL)
## dba.plotProfile(profiles)


###################################################
### code chunk number 41: DiffBind.Rnw:1083-1084
###################################################
tamoxifen <- dba.contrast(tamoxifen,design="~Tissue + Condition")


###################################################
### code chunk number 42: DiffBind.Rnw:1090-1092
###################################################
tamoxifen <- dba.analyze(tamoxifen)
dba.show(tamoxifen, bContrasts=TRUE)


###################################################
### code chunk number 43: tamox_block_ma
###################################################
dba.plotMA(tamoxifen)


###################################################
### code chunk number 44: tamox_block_vol
###################################################
dba.plotVolcano(tamoxifen)


###################################################
### code chunk number 45: DiffBind.Rnw:1141-1142
###################################################
multifactor.DB <- dba.report(tamoxifen)


###################################################
### code chunk number 46: DiffBind.Rnw:1149-1151
###################################################
min(abs(tamoxifen.DB$Fold))
min(abs(multifactor.DB$Fold))


###################################################
### code chunk number 47: DiffBind.Rnw:1159-1161
###################################################
sum(tamoxifen.DB$Fold > 0) / sum(tamoxifen.DB$Fold < 0)
sum(multifactor.DB$Fold > 0) / sum(multifactor.DB$Fold < 0)


###################################################
### code chunk number 48: DiffBind.Rnw:1167-1169
###################################################
tamoxifen <- dba.analyze(tamoxifen,method=DBA_ALL_METHODS)
dba.show(tamoxifen,bContrasts=TRUE)


###################################################
### code chunk number 49: tamox_block_venn
###################################################
tamoxifen.OL <- dba.plotVenn(tamoxifen,contrast=1,method=DBA_ALL_METHODS,
                             bDB=TRUE)


###################################################
### code chunk number 50: DiffBind.Rnw:1209-1213
###################################################
tamoxifen <- dba.contrast(tamoxifen,contrast=c("Tissue","MCF7","T47D"),
                          reorderMeta = list(Tissue="MCF7"))
tamoxifen <- dba.analyze(tamoxifen,method=DBA_ALL_METHODS)
dba.show(tamoxifen,bContrasts=TRUE)


###################################################
### code chunk number 51: blacklistPeaks
###################################################
data(tamoxifen_peaks)
tamoxifen
peakdata  <- dba.show(tamoxifen)$Intervals
tamoxifen <- dba.blacklist(tamoxifen, blacklist=DBA_BLACKLIST_HG19, 
                           greylist=FALSE)
tamoxifen
peakdata.BL <- dba.show(tamoxifen)$Intervals
peakdata - peakdata.BL


###################################################
### code chunk number 52: blacklistAnal
###################################################
length(multifactor.DB)
data(tamoxifen_counts)
tamoxifen   <- dba.blacklist(tamoxifen, blacklist=DBA_BLACKLIST_HG19,
                             greylist=FALSE)
blacklisted <- dba.blacklist(tamoxifen, Retrieve=DBA_BLACKLISTED_PEAKS)
tamoxifen   <- dba.contrast(tamoxifen, design="~Tissue + Condition")
tamoxifen   <- dba.analyze(tamoxifen)
blacklisted.DB <- dba.report(tamoxifen)
length(blacklisted.DB)


###################################################
### code chunk number 53: blacklistRes
###################################################
bl_site <- match(blacklisted[[1]], multifactor.DB)
multifactor.DB[bl_site,]
is.na(match(blacklisted[[1]], blacklisted.DB))


###################################################
### code chunk number 54: greylistGet
###################################################
data(tamoxifen_greylist)
names(tamoxifen.greylist)
tamoxifen.greylist$master


###################################################
### code chunk number 55: greylistControls
###################################################
names(tamoxifen.greylist$controls)
tamoxifen.greylist$controls


###################################################
### code chunk number 56: greylistPeaks
###################################################
data(tamoxifen_peaks)
tamoxifen <- dba.blacklist(tamoxifen, blacklist=DBA_BLACKLIST_HG19,
                           greylist=tamoxifen.greylist)


###################################################
### code chunk number 57: greylistCons
###################################################
data(tamoxifen_counts)
cons.peaks <- dba.show(tamoxifen)$Intervals[1]
tamoxifen  <- dba.blacklist(tamoxifen, blacklist=DBA_BLACKLIST_HG19,
                            greylist=tamoxifen.greylist)
cons.peaks.grey <- dba.show(tamoxifen)$Intervals[1]
cons.peaks - cons.peaks.grey


###################################################
### code chunk number 58: greylistMake (eval = FALSE)
###################################################
## tamoxifen <- dba(sampleSheet="tamoxifen.csv")
## tamoxifen <- dba.blacklist(tamoxifen)
## tamoxifen.greylist <- dba.blacklist(tamoxifen, Retrieve=DBA_GREYLIST)


###################################################
### code chunk number 59: greylistP
###################################################
tamoxifen$config$greylist.pval <- 0.999


###################################################
### code chunk number 60: norm0
###################################################
data(tamoxifen_analysis)
dba.plotMA(tamoxifen, contrast=list(Resistant=tamoxifen$masks$Resistant),
           bNormalized=FALSE, sub="Non-Normalized")


###################################################
### code chunk number 61: normDESeq2LibFull
###################################################
tamoxifen <- dba.normalize(tamoxifen, normalize=DBA_NORM_LIB)
tamoxifen <- dba.analyze(tamoxifen)
dba.plotMA(tamoxifen, method=DBA_DESEQ2, sub="DESeq2:lib:full")


###################################################
### code chunk number 62: normRes1
###################################################
dbs <- dba.report(tamoxifen, bDB=TRUE, bGain=TRUE, bLoss=TRUE)
dbs$config$factor <- "normalize"
dbs$class[DBA_ID,] <- colnames(dbs$class)[1] <-  "LIB_Full"
dbs$class[DBA_FACTOR,] <- DBA_NORM_LIB
dbs


###################################################
### code chunk number 63: normDESeq2RLE
###################################################
tamoxifen <- dba.normalize(tamoxifen, normalize=DBA_NORM_NATIVE)
tamoxifen <- dba.analyze(tamoxifen)
dba.plotMA(tamoxifen, method=DBA_DESEQ2, sub="DESeq2:RLE:RiP")


###################################################
### code chunk number 64: normRes2
###################################################
db <- dba.report(tamoxifen, bDB=TRUE, bGain=TRUE, bLoss=TRUE)
db$class[DBA_ID,] <- "RLE_RiP"
db$class[DBA_FACTOR,] <- DBA_NORM_RLE
dbs <- dba.peakset(dbs, db)
db


###################################################
### code chunk number 65: normDESeq2Comparison
###################################################
par(mfrow=c(3,1))
dba.plotVenn(dbs,c(1,4), main="Total DB Sites")
dba.plotVenn(dbs,dbs$masks$Gain,main="Gain in Resistant")
dba.plotVenn(dbs,dbs$masks$Loss,main="Gain in Responsive")
par(mfrow=c(1,1))


###################################################
### code chunk number 66: normDESeq2LibRiP
###################################################
tamoxifen <- dba.normalize(tamoxifen, normalize=DBA_NORM_LIB,
                           library=DBA_LIBSIZE_PEAKREADS, background=FALSE)
tamoxifen <- dba.analyze(tamoxifen)
dba.plotMA(tamoxifen, method=DBA_DESEQ2, sub="DESeq2:lib:RiP")


###################################################
### code chunk number 67: normRes3
###################################################
dbs$class[DBA_CONDITION,1:3] <- DBA_LIBSIZE_FULL
dbs$class[DBA_CONDITION,4:6] <- DBA_LIBSIZE_PEAKREADS
dbs$config$condition <- "lib.size"
db <- dba.report(tamoxifen, bDB=TRUE, bGain=TRUE, bLoss=TRUE)
db$class[DBA_ID,] <- "LIB_RiP"
db$class[DBA_FACTOR,] <- DBA_NORM_LIB
db$class[DBA_CONDITION,] <- DBA_LIBSIZE_PEAKREADS
dbs <- dba.peakset(dbs, db)
db


###################################################
### code chunk number 68: normDESeq2LibsizeVenn
###################################################
dba.plotVenn(dbs,c(1,7,4),main="DB Sites")


###################################################
### code chunk number 69: frip
###################################################
dba.show(tamoxifen,attributes=c(DBA_ID,DBA_FRIP))


###################################################
### code chunk number 70: normBG
###################################################
data(tamoxifen_analysis)
tamoxifen <- dba.normalize(tamoxifen, method=DBA_ALL_METHODS,
                           normalize=DBA_NORM_NATIVE,
                           background=TRUE)
tamoxifen <- dba.analyze(tamoxifen, method=DBA_ALL_METHODS)
dba.show(tamoxifen,bContrasts=TRUE)
par(mfrow=c(2,1))
dba.plotMA(tamoxifen, method=DBA_EDGER, sub="edgeR:TMM:background")
dba.plotMA(tamoxifen, method=DBA_DESEQ2, sub="DESeq2:RLE:background")
par(mfrow=c(1,1))


###################################################
### code chunk number 71: normBgRes
###################################################
db <- dba.report(tamoxifen, bDB=TRUE, bGain=TRUE, bLoss=TRUE)
db$class[DBA_ID,] <- "RLE_BG"
db$class[DBA_FACTOR,] <- DBA_NORM_RLE
db$class[DBA_CONDITION,] <- DBA_LIBSIZE_BACKGROUND
dbs <- dba.peakset(dbs, db)
db


###################################################
### code chunk number 72: normBGVenns
###################################################
par(mfcol=c(3,2))
dba.plotVenn(dbs,c(1,10),   main="All Differentially Bound Sites")
dba.plotVenn(dbs,c(2,11),   main="Gain in Resistant cells")
dba.plotVenn(dbs,c(3,12),   main="Loss in Resistant cells")
dba.plotVenn(dbs,c(1,10,4), main="All Differentially Bound Sies")
dba.plotVenn(dbs,c(2,11,5), main="Gain in Resistant cells")
dba.plotVenn(dbs,c(3,12,6), main="Loss in Resistant cells")


###################################################
### code chunk number 73: MCF7T47D
###################################################
mcf7t47d <- dba(tamoxifen,mask=c(3:7))
dba.plotMA(mcf7t47d,
           contrast=list(MCF7=mcf7t47d$masks$MCF7,
                         T47D=mcf7t47d$masks$T47D), 
           bNormalized=FALSE)


###################################################
### code chunk number 74: loessfit
###################################################
mcf7t47d$config$AnalysisMethod <- DBA_EDGER
mcf7t47d <- dba.normalize(mcf7t47d, offsets=TRUE)
mcf7t47d <- dba.contrast(mcf7t47d, contrast=c("Tissue","MCF7","T47D"))
mcf7t47d <- dba.analyze(mcf7t47d)
dba.plotMA(mcf7t47d)


###################################################
### code chunk number 75: compareLoess
###################################################
mcf7t47d.DB <- dba.report(mcf7t47d)
sum(mcf7t47d.DB$Fold > 0)
sum(mcf7t47d.DB$Fold < 0)


###################################################
### code chunk number 76: normCompare
###################################################
data(tamoxifen_analysis)
dbs.all <- NULL
for(norm in c("lib","RLE","TMM", "loess")) {
  for(libsize in c("full","RiP","background")) {
    tam <- NULL
    background <- offsets <- FALSE
    if(libsize == "full" && norm != "lib") {
      background <- NULL
    }
    if(libsize == DBA_LIBSIZE_BACKGROUND) {
      background <- TRUE
      if(norm == DBA_NORM_LIB) {
        background <- NULL
      }
    } 
    if(norm == "loess" && !is.null(background)) {
      offsets <- TRUE
      if(libsize != "background") {
        background <- FALSE
      } else {
        background <- NULL
      }
    }
    if(!is.null(background)) {
      tam <- dba.normalize(tamoxifen, method=DBA_ALL_METHODS,
                           normalize=norm, library=libsize, 
                           background=background, offsets=offsets)
    }
    if(!is.null(tam)) {
      tam <- dba.analyze(tam, method=DBA_ALL_METHODS)
      for(meth in DBA_ALL_METHODS) {
        db <- dba.report(tam, method=meth, bDB=TRUE)
        if(meth == DBA_EDGER) {
          methstr <- "edgeR"
        } else {
          methstr <- "DESeq2"
        }
        if(libsize == "background") {
          libstr <- "BG"
        } else {
          libstr <- libsize
        }
        id <- paste(norm,libstr,methstr,sep="_")
        if(libsize == "full") {
          libstr <- "BG"
        }
        if(is.null(dbs.all)) {
          dbs.all <- db
          dbs.all$config$factor    <- "Normalization Method"
          dbs.all$config$condition <- "Reference Reads"
          dbs.all$config$treatment <- "Analysis Method"
          dbs.all$class[DBA_ID,]  <- colnames(dbs.all$class)[1] <- id
          
          dbs.all$class[DBA_FACTOR,]    <- norm
          dbs.all$class[DBA_CONDITION,] <- libstr
          dbs.all$class[DBA_TREATMENT,] <- "edgeR"
          dbs.all$class[DBA_TISSUE,]    <- NA
        } else {
          db$class[DBA_ID,]        <- id
          db$class[DBA_FACTOR,]    <- norm
          db$class[DBA_CONDITION,] <- libstr
          db$class[DBA_TREATMENT,] <- methstr
          db$class[DBA_TISSUE,]    <- NA
          dbs.all <- dba.peakset(dbs.all,db)
        }
      }
    }
  }
}
dbs.all <- dba(dbs.all,minOverlap=1)
dbs.all


###################################################
### code chunk number 77: normClusterDESeq2
###################################################
deseq <- dba(dbs.all,mask=dbs.all$masks$DESeq2, minOverlap = 1)
binding <- dba.peakset(deseq, bRetrieve=TRUE)
dba.plotHeatmap(deseq, maxSites=nrow(binding),  bLog=FALSE,
                correlations=FALSE,minval=-5, maxval=5, cexCol=1.3, 
                colScheme = hmap, main = "DESeq2 Differentially Bound Sites", 
                ColAttributes = c(DBA_CONDITION, DBA_FACTOR),
                key.title = "LFC")


###################################################
### code chunk number 78: normCluster1
###################################################
binding <- dba.peakset(dbs.all, bRetrieve=TRUE)
dba.plotHeatmap(dbs.all, maxSites=nrow(binding), bLog=FALSE,
                correlations=FALSE, minval=-5, maxval=5, cexCol=1.3, 
                colScheme = hmap, key.title="LFC",
                ColAttributes = c(DBA_CONDITION, DBA_TREATMENT, DBA_FACTOR),
                main="All Differentially Bound Sites")


###################################################
### code chunk number 79: normCluster2
###################################################
dba.plotHeatmap(dbs.all, cexCol=1.3, main="Correlations of DB Sites",
                ColAttributes = c(DBA_CONDITION, DBA_TREATMENT, DBA_FACTOR))


###################################################
### code chunk number 80: loadSpikes
###################################################
load(system.file('extra/spikes.rda',package='DiffBind'))
spikes


###################################################
### code chunk number 81: spikeins
###################################################
spikes$samples$Spikein


###################################################
### code chunk number 82: spikeMAnone
###################################################
dba.plotMA(spikes, contrast=list(Fulvestrant=spikes$masks$Fulvestrant),
           bNormalized=FALSE, sub="RAW", bSmooth=FALSE, dotSize=1.5)


###################################################
### code chunk number 83: spikeMAstraight
###################################################
par(mfrow=c(3,1))
spikes <- dba.normalize(spikes, normalize=DBA_NORM_LIB, 
                        background=FALSE)
spikes <- dba.analyze(spikes)
dba.plotMA(spikes, sub="LIB full", bSmooth=FALSE, dotSize=1.5)

spikes <- dba.normalize(spikes, normalize="RLE",
                        background=FALSE)
spikes <- dba.analyze(spikes)
dba.plotMA(spikes, sub="RLE RiP", bSmooth=FALSE, dotSize=1.5)

spikes <- dba.normalize(spikes, normalize="RLE", 
                        background=TRUE)
spikes <- dba.analyze(spikes)
dba.plotMA(spikes, sub="RLE BG", bSmooth=FALSE, dotSize=1.5)
par(mfrow=c(1,1))


###################################################
### code chunk number 84: spikeMAspikein
###################################################
par(mfrow=c(2,1))
spikes <- dba.normalize(spikes, normalize=DBA_NORM_LIB, 
                        spikein=spikes.spikeins)
spikes <- dba.analyze(spikes)
dba.plotMA(spikes, sub="LIB spikein", bSmooth=FALSE, dotSize=1.5)

spikes <- dba.normalize(spikes, normalize=DBA_NORM_RLE, spikein = TRUE)
spikes <- dba.analyze(spikes)
dba.plotMA(spikes, sub="RLE spikein", bSmooth=FALSE, dotSize=1.5)
par(mfrow=c(1,1))


###################################################
### code chunk number 85: loadParfac
###################################################
load(system.file('extra/parallelFactor.rda',package='DiffBind'))
parallelFactor


###################################################
### code chunk number 86: parfacMAnone
###################################################
dba.plotMA(parallelFactor,
           contrast=list(Fulvestrant=parallelFactor$masks$Fulvestrant),
           bNormalized=FALSE, sub="RAW", bSmooth=FALSE, dotSize=1.5)


###################################################
### code chunk number 87: parfacMA
###################################################
par(mfrow=c(2,1))
parallelFactor <- dba.normalize(parallelFactor, norm=DBA_NORM_LIB,
                                spikein = parallelFactor.peaks)
parallelFactor <- dba.analyze(parallelFactor)
dba.plotMA(parallelFactor, sub="LIB CTCF", bSmooth=FALSE, dotSize=1.5)

parallelFactor <- dba.normalize(parallelFactor, norm=DBA_NORM_RLE,
                                spikein = TRUE) 
parallelFactor <- dba.analyze(parallelFactor)
dba.plotMA(parallelFactor, sub="RLE CTCF", bSmooth=FALSE, dotSize=1.5)
par(mfrow=c(1,1))


###################################################
### code chunk number 88: DiffBind.Rnw:2385-2411
###################################################
require(xtable)
collabs <- paste(DBA_NORM_LIB," & ",DBA_NORM_RLE," & ",
                 DBA_NORM_TMM," & ",DBA_OFFSETS_LOESS)

rowlabs <- c(DBA_LIBSIZE_FULL,DBA_LIBSIZE_PEAKREADS,
             DBA_LIBSIZE_BACKGROUND,
             DBA_NORM_SPIKEIN,"parallel factor")

normtab <- matrix("X",5,4)
normtab[1,2:4] <- ""
normtab[c(1,3:5),4] <- ""
normtab <- data.frame(normtab)
rownames(normtab) <- rowlabs
addtorow <- list()
addtorow$pos <- list(0,0)
addtorow$command <- c("& \\multicolumn{4}{c|}{Normalization} \\\\\n",
                      paste("Reference & ",collabs,"\\\\\n"))

captionStr <- paste("Table of allowable normalization schemes.",
                    "Columns are normalization methods set by \\Rcode{normalize} or \\Rcode{offsets}.",
                    "Rows are reference reads, set by \\Rcode{library}, \\Rcode{background}, or \\Rcode{spikein}.")

print(xtable::xtable(normtab,align=c(rep("|c",5),"|"),
                     caption=captionStr, label="DiffBind-normtable"), 
      add.to.row = addtorow,
      include.colnames=FALSE, scalebox=1)


###################################################
### code chunk number 89: DiffBind.Rnw:2452-2453
###################################################
data(tamoxifen_peaks)


###################################################
### code chunk number 90: DiffBind.Rnw:2477-2479
###################################################
olap.rate <- dba.overlap(tamoxifen,mode=DBA_OLAP_RATE)
olap.rate


###################################################
### code chunk number 91: tamox_rate
###################################################
plot(olap.rate,type='b',ylab='# peaks',
     xlab='Overlap at least this many peaksets')


###################################################
### code chunk number 92: DiffBind.Rnw:2536-2537
###################################################
names(tamoxifen$masks)


###################################################
### code chunk number 93: DiffBind.Rnw:2546-2548
###################################################
dba.overlap(tamoxifen,tamoxifen$masks$MCF7 & tamoxifen$masks$Responsive,
            mode=DBA_OLAP_RATE)


###################################################
### code chunk number 94: tamox_mcf7_venn
###################################################
dba.plotVenn(tamoxifen, tamoxifen$masks$MCF7 & tamoxifen$masks$Responsive)


###################################################
### code chunk number 95: DiffBind.Rnw:2579-2582
###################################################
tamoxifen_consensus <- dba.peakset(tamoxifen, 
                                   consensus=c(DBA_TISSUE,DBA_CONDITION),
                                   minOverlap=0.66)


###################################################
### code chunk number 96: DiffBind.Rnw:2597-2601
###################################################
tamoxifen_consensus <- dba(tamoxifen_consensus,
                           mask=tamoxifen_consensus$masks$Consensus,
                           minOverlap=1)
tamoxifen_consensus


###################################################
### code chunk number 97: DiffBind.Rnw:2607-2608
###################################################
consensus_peaks <- dba.peakset(tamoxifen_consensus, bRetrieve=TRUE)


###################################################
### code chunk number 98: tamox_lines_venn
###################################################
data(tamoxifen_peaks)
tamoxifen <- dba.peakset(tamoxifen, consensus=DBA_TISSUE, minOverlap=0.66)
cons.ol <- dba.plotVenn(tamoxifen, tamoxifen$masks$Consensus)


###################################################
### code chunk number 99: DiffBind.Rnw:2648-2649
###################################################
data(tamoxifen_peaks)


###################################################
### code chunk number 100: DiffBind.Rnw:2655-2657
###################################################
dba.overlap(tamoxifen,tamoxifen$masks$Resistant,mode=DBA_OLAP_RATE)
dba.overlap(tamoxifen,tamoxifen$masks$Responsive,mode=DBA_OLAP_RATE)


###################################################
### code chunk number 101: tamox_cons_venn
###################################################
tamoxifen <- dba.peakset(tamoxifen, consensus=DBA_CONDITION, minOverlap=0.33)
dba.plotVenn(tamoxifen,tamoxifen$masks$Consensus)


###################################################
### code chunk number 102: DiffBind.Rnw:2692-2693
###################################################
tamoxifen.OL <- dba.overlap(tamoxifen, tamoxifen$masks$Consensus)


###################################################
### code chunk number 103: DiffBind.Rnw:2700-2702
###################################################
tamoxifen.OL$onlyA
tamoxifen.OL$onlyB


###################################################
### code chunk number 104: tamox_compare_venn
###################################################
tamoxifen <- dba.peakset(tamoxifen,tamoxifen$masks$Consensus,
                         minOverlap=1,sampID="OL Consensus")
tamoxifen <- dba.peakset(tamoxifen,!tamoxifen$masks$Consensus,
                         minOverlap=3,sampID="Consensus_3")
dba.plotVenn(tamoxifen,14:15)


###################################################
### code chunk number 105: DiffBind.Rnw:2738-2739
###################################################
data(tamoxifen_analysis)


###################################################
### code chunk number 106: DiffBind.Rnw:2749-2750
###################################################
tamoxifen.rep <- dba.report(tamoxifen,bCalled=TRUE,th=1)


###################################################
### code chunk number 107: DiffBind.Rnw:2759-2765
###################################################
onlyResistant <- tamoxifen.rep$Called1>=2 & tamoxifen.rep$Called2<3
sum(onlyResistant )
onlyResponsive <- tamoxifen.rep$Called2>=3 &  tamoxifen.rep$Called1<2
sum(onlyResponsive)
bothGroups <- tamoxifen.rep$Called1>= 2 & tamoxifen.rep$Called2>=3
sum(bothGroups)


###################################################
### code chunk number 108: DiffBind.Rnw:2780-2789
###################################################
tamoxifen.DB <- dba.report(tamoxifen,bCalled=TRUE)
onlyResistant.DB  <- (tamoxifen.DB$Called1 >= 2) & (tamoxifen.DB$Called2 < 3)
sum(onlyResistant.DB)
onlyResponsive.DB <- (tamoxifen.DB$Called2 >= 3) & (tamoxifen.DB$Called1 < 2)
sum(onlyResponsive.DB)
bothGroups.DB     <- (tamoxifen.DB$Called1 >= 2) & (tamoxifen.DB$Called2 >= 3)
sum(bothGroups.DB)
neitherGroup.DB   <- (tamoxifen.DB$Called1  < 2) & (tamoxifen.DB$Called2 < 3)
sum(neitherGroup.DB)


###################################################
### code chunk number 109: DiffBind.Rnw:2892-2893 (eval = FALSE)
###################################################
## DBA$config$design <- FALSE


###################################################
### code chunk number 110: DiffBind.Rnw:3307-3314 (eval = FALSE)
###################################################
## tmpdir <- tempdir()
## url <- 'https://content.cruk.cam.ac.uk/bioinformatics/software/DiffBind/DiffBind_vignette_data.tar.gz'
## file <- basename(url)
## options(timeout=600)
## download.file(url, file.path(tmpdir,file))
## untar(file.path(tmpdir,file), exdir = tmpdir )
## setwd(file.path(tmpdir,"DiffBind_Vignette"))


###################################################
### code chunk number 111: sessionInfo
###################################################
toLatex(sessionInfo())


