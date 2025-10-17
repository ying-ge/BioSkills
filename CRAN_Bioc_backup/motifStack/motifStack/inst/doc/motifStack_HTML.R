## ----checkplatform, include=FALSE---------------------------------------------
library(knitr)
library(motifStack)
library(MotifDb)
library(ade4)
library(ggplot2)
library(TFBSTools)
library(JASPAR2020)
if(.Platform$OS.type=="windows"){
  opts_chunk$set(eval=FALSE)
}

## ----importMatrix,fig.cap="import data from PFMatrixList",fig.width=6,fig.height=3----
library(motifStack)
library(JASPAR2020)
motifs <- importMatrix(getMatrixSet(JASPAR2020, 
                                    list(species="Mus musculus")))
plot(motifs[[1]])

## ----importFromFiles,fig.cap="import data from file",fig.width=6,fig.height=3----
RUNX1 <- importMatrix(system.file("extdata", "MA0002.1.jaspar",
                                  package = "motifStack",
                                  mustWork = TRUE))[[1]]
plot(RUNX1)

## ----DNAseqLogo,fig.cap="Plot a DNA sequence logo with different fonts and colors",fig.width=8,fig.height=2.5----
library(motifStack)
pcm <- read.table(file.path(find.package("motifStack"), 
                            "extdata", "bin_SOLEXA.pcm"))
pcm <- pcm[,3:ncol(pcm)]
rownames(pcm) <- c("A","C","G","T")
motif <- new("pcm", mat=as.matrix(pcm), name="bin_SOLEXA")
##pfm object
#motif <- pcm2pfm(pcm)
#motif <- new("pfm", mat=motif, name="bin_SOLEXA")
plot(motif)
#plot the logo with same height
plot(motif, ic.scale=FALSE, ylab="probability")
#try a different font and a different color group
motif@color <- colorset(colorScheme='basepairing')
plot(motif,font="serif")

## ----seqLogoMarkers, fig.cap="Plot a DNA sequence logo with markers", fig.width=8,fig.height=2.5----
markerRect <- new("marker", type="rect", start=6, stop=7, gp=gpar(lty=2, fill=NA, col="orange"))
markerLine <- new("marker", type="line", start=2, stop=7, gp=gpar(lwd=2, col="red"))
markerText <- new("marker", type="text", start=c(1, 5), 
                  label=c("*", "core"), gp=gpar(cex=2, col="red"))
motif <- new("pcm", mat=as.matrix(pcm), name="bin_SOLEXA", 
             markers=c(markerRect, markerLine, markerText))
plot(motif)

## ----seqLogoXaxis, fig.cap="Plot a DNA sequence logo with pre-defined xlabels", fig.width=8,fig.height=2.5----
plot(motif, xaxis=paste0("pos", seq.int(7)+10))

## ----RNAseqLogo,fig.cap="Plot an RNA sequence logo",fig.width=6,fig.height=3----
rna <- pcm
rownames(rna)[4] <- "U"
motif <- new("pcm", mat=as.matrix(rna), name="RNA_motif")
plot(motif)

## ----AAseqLogo,fig.cap="Plot an sequence logo with any symbols as you want such as amino acid sequence logo",fig.width=6,fig.height=3----
library(motifStack)
protein<-read.table(file.path(find.package("motifStack"),"extdata","cap.txt"))
protein<-t(protein[,1:20])
motif<-pcm2pfm(protein)
motif<-new("pfm", mat=motif, name="CAP", 
            color=colorset(alphabet="AA",colorScheme="chemistry"))
plot(motif)

## ----affinityLogo,fig.cap="Plot an affinity logo",fig.width=6,fig.height=3----
library(motifStack)
motif<-matrix(
  c(
    .846, .631, .593, .000, .000, .000, .434, .410, 1.00, .655, .284, .000, .000, .771, .640, .961,
    .625, .679, .773, 1.00, 1.00, .000, .573, .238, .397, 1.00, 1.00, .000, .298, 1.00, 1.00, .996,
    1.00, 1.00, 1.00, .228, .000, 1.00, 1.00, .597, .622, .630, .000, 1.00, 1.00, .871, .617, 1.00,
    .701, .513, .658, .000, .000, .247, .542, 1.00, .718, .686, .000, .000, .000, .595, .437, .970
  ), nrow=4, byrow = TRUE)
rownames(motif) <- c("A", "C", "G", "T")
motif<-new("psam", mat=motif, name="affinity logo", 
           markers=list(new("marker", type="rect",
                            start=c(4, 11), stop=c(6, 13),
                            gp=gpar(col="#009E73", fill=NA, lty=2))))
plot(motif)

## ----logostack,fig.cap="Plot motifs with sequence logo stack style",fig.width=4,fig.height=6----
library(motifStack)
#####Input#####
motifs<-importMatrix(dir(file.path(find.package("motifStack"),
                                   "extdata"),"pcm$", 
                         full.names = TRUE))

## plot stacks
motifStack(motifs, layout="stack", ncex=1.0)

## ----rnalogostack,fig.cap="Plot RNA motifs with sequence logo stack style", fig.width=4,fig.height=6----
rnaMotifs <- DNAmotifToRNAmotif(motifs)
names(rnaMotifs)
motifStack(rnaMotifs, layout = "stack", 
           reorder=FALSE) ## we can also use reorder=FALSE to keep the order of input. 

## ----logostack2,fig.cap="Plot affinity logos with sequence logo stack style",fig.width=4,fig.height=3.5----
motif2 <- motif
motif2$mat <- motif$mat[, 5:12]
motif2$name <- "logo2"
psamMotifs <- list(motif, motif2)
motifStack(psamMotifs)

## ----treestack,fig.cap="Sequence logo stack with hierarchical cluster tree",fig.width=5,fig.height=6----
## plot stacks with hierarchical tree
motifStack(motifs, layout="tree")

## ----radialstack,fig.cap="Plot motifs in a radial style when the number of motifs is too much to be shown in a vertical stack",fig.width=6,fig.height=6----
## When the number of motifs is too much to be shown in a vertical stack, 
## motifStack can draw them in a radial style.
## random sample from MotifDb
library("MotifDb")
matrix.fly <- query(MotifDb, "Dmelanogaster")
motifs2 <- as.list(matrix.fly)
## use data from FlyFactorSurvey
motifs2 <- motifs2[grepl("Dmelanogaster\\-FlyFactorSurvey\\-",
                         names(motifs2))]
## format the names
names(motifs2) <- gsub("Dmelanogaster_FlyFactorSurvey_", "",
                       gsub("_FBgn\\d+$", "",
                            gsub("[^a-zA-Z0-9]","_",
                                 gsub("(_\\d+)+$", "", names(motifs2)))))
motifs2 <- motifs2[unique(names(motifs2))]
pfms <- sample(motifs2, 30)
## creat a list of object of pfm 
motifs2 <- mapply(pfms, names(pfms), FUN=function(.ele, .name){
  new("pfm",mat=.ele, name=.name)}, SIMPLIFY = FALSE)
## trim the motifs
motifs2 <- lapply(motifs2, trimMotif, t=0.4)
## setting colors
library(RColorBrewer)
color <- brewer.pal(10, "Set3")
## plot logo stack with radial style
motifStack(motifs2, layout="radialPhylog", 
           circle=0.3, cleaves = 0.2, 
           clabel.leaves = 0.5, 
           col.bg=rep(color, each=3), col.bg.alpha=0.3, 
           col.leaves=rep(color, each=3),
           col.inner.label.circle=rep(color, each=3), 
           inner.label.circle.width=0.05,
           col.outer.label.circle=rep(color, each=3), 
           outer.label.circle.width=0.02, 
           circle.motif=1.2,
           angle=350)

## ----motifCloud,fig.cap="Sequence logo cloud with rectangle packing layout",fig.width=6,fig.height=6----
## assign groups for motifs
groups <- rep(paste("group",1:5,sep=""), each=10)
names(groups) <- names(pfms)
## assign group colors
group.col <- brewer.pal(5, "Set3")
names(group.col)<-paste("group",1:5,sep="")
## create a list of pfm objects
pfms <- mapply(names(pfms), pfms, FUN=function(.ele, .pfm){
  new("pfm",mat=.pfm, name=.ele)}
               ,SIMPLIFY = FALSE)
## use matalign to calculate the distances of motifs
hc <- clusterMotifs(pfms)
## convert the hclust to phylog object
library(ade4)
phylog <- ade4::hclust2phylog(hc)
## reorder the pfms by the order of hclust
leaves <- names(phylog$leaves)
pfms <- pfms[leaves]
## extract the motif signatures
motifSig <- motifSignature(pfms, phylog, cutoffPval=0.0001, min.freq=1)
## draw the motifs with a tag-cloud style.
motifCloud(motifSig, scale=c(6, .5), 
           layout="rectangles", 
           group.col=group.col, 
           groups=groups, 
           draw.legend=TRUE)

## ----motifRadialPhylog,fig.cap="Grouped sequence logo with radialPhylog style layout",fig.width=6,fig.height=6----
## get the signatures from object of motifSignature
sig <- signatures(motifSig)
## set the inner-circle color for each signature
gpCol <- sigColor(motifSig)
## plot the logo stack with radial style.
plotMotifStackWithRadialPhylog(phylog=phylog, pfms=sig, 
                              circle=0.4, cleaves = 0.3, 
                              clabel.leaves = 0.5, 
                              col.bg=rep(color, each=3), col.bg.alpha=0.3, 
                              col.leaves=rep(rev(color), each=3),
                              col.inner.label.circle=gpCol, 
                              inner.label.circle.width=0.03,
                              angle=350, circle.motif=1.2, 
                              motifScale="logarithmic")

## ----motifCircos,fig.cap="Grouped sequence logo with circos style layout",fig.width=6,fig.height=6----
## plot the logo stack with cirsoc style.
motifCircos(phylog=phylog, pfms=pfms, pfms2=sig, 
            col.tree.bg=rep(color, each=5), col.tree.bg.alpha=0.3, 
            col.leaves=rep(rev(color), each=5),
            col.inner.label.circle=gpCol, 
            inner.label.circle.width=0.03,
            col.outer.label.circle=gpCol, 
            outer.label.circle.width=0.03,
            r.rings=c(0.02, 0.03, 0.04), 
            col.rings=list(sample(colors(), 30), 
                           sample(colors(), 30), 
                           sample(colors(), 30)),
            angle=350, motifScale="logarithmic")

## ----motifPilesHeatmap,fig.cap="Grouped sequence logo with a heatmap",fig.width=6,fig.height=6----
## plot the logo stack with heatmap.
df <- data.frame(A=runif(n = 30), B=runif(n = 30), C=runif(n = 30), D=runif(n = 30))
map2col <- function(x, pal){
  rg <- range(x)
  pal[findInterval(x, seq(rg[1], rg[2], length.out = length(pal)+1), 
                   all.inside = TRUE)]
}
dl <- lapply(df, map2col, pal=heat.colors(10))
## alignment of the pfms, this step will make the motif logos occupy 
## more space. Users can skip this alignment to see the difference.
pfmsAligned <- DNAmotifAlignment(pfms)
## plot motifs
motifPiles(phylog=phylog, pfms=pfmsAligned, 
            col.tree=rep(color, each=5),
            col.leaves=rep(rev(color), each=5),
            col.pfms2=gpCol, 
            r.anno=rep(0.02, length(dl)), 
            col.anno=dl,
            motifScale="logarithmic",
            plotIndex=TRUE,
            groupDistance=10)

## ----browseMotifs,fig.width=8,fig.height=8------------------------------------
browseMotifs(pfms = pfms, phylog = phylog, layout="tree", yaxis = FALSE, baseWidth=6, baseHeight = 15)

## ----browseMotifsRadialPhylog,fig.width=8,fig.height=8------------------------
browseMotifs(pfms = pfms, phylog = phylog, layout="radialPhylog", yaxis = FALSE, xaxis = FALSE, baseWidth=6, baseHeight = 15)

## ----geommotif,fig.width=8, fig.height=8--------------------------------------
pcm <- read.table(file.path(find.package("motifStack"), 
                            "extdata", "bin_SOLEXA.pcm"))
pcm <- pcm[,3:ncol(pcm)]
rownames(pcm) <- c("A","C","G","T")
markerRect <- new("marker", type="rect", start=6, stop=7, gp=gpar(lty=2, fill=NA, col="orange"))
markerLine <- new("marker", type="line", start=3, stop=5, gp=gpar(lwd=2, col="red"))
markerText <- new("marker", type="text", start=1, label="*", gp=gpar(cex=2, col="red"))
motif <- new("pcm", mat=as.matrix(pcm), name="bin_SOLEXA", 
             markers=c(markerRect, markerLine, markerText))
pfm <- pcm2pfm(motif)
df <- data.frame(xmin=c(.25, .25), ymin=c(.25, .75), xmax=c(.75, .75), ymax=c(.5, 1), 
                 fontfamily=c("serif", "mono"), fontface=c(2, 1))
df$motif <- list(pfm, pfm)

library(ggplot2)

ggplot(df, aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, motif=motif, 
               fontfamily=fontfamily, fontface=fontface)) + 
    geom_motif() + theme_bw() + ylim(0, 1) + xlim(0, 1)

df <- data.frame(x=.5, y=c(.25, .75), width=.5, height=.25, 
                 fontfamily=c("serif", "mono"), fontface=c(2, 1))
df$motif <- list(pfm, pfm)

ggplot(df, aes(x=x, y=y, width=width, height=height, motif=motif, 
               fontfamily=fontfamily, fontface=fontface)) + 
    geom_motif(use.xy=TRUE) + theme_bw() + ylim(0, 1) + xlim(0, 1)

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

