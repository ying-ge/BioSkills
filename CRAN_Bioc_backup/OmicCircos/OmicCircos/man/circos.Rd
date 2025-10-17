%% circos
\name{circos}
\alias{circos}

\title{
 draw circular
}
\description{
 This is the main function of OmicCircos to draw circular plots. 
}
\usage{
  circos(mapping=mapping, xc=400, yc=400, R=400, W=W,  cir="", type="n", 
    col.v=3, B=FALSE, print.chr.lab=TRUE, col.bar=FALSE, col.bar.po = "topleft",       
    cluster=FALSE, order=NULL, scale=FALSE, cutoff = "n", zoom="", cex=1, lwd=1, 
    col=rainbow(10, alpha=0.8)[7], side="")
}

\arguments{
	\item{mapping}{data frame or matrix containing mapping information 
  and values. Column 1: the segment or chromosome ID; column 2: the position; 
  column 3: the position label (optional) or the value and additional columns are the values.   
such as gene expression and copy number.  Missing values are allowed and will be ignored.}
  
   \item{xc}{integer, the circle center x coordinate}
   \item{yc}{integer, the circle center y coordinate}
   \item{R}{integer, the circle radius}
   \item{W}{integer, the circle width}
   \item{cir}{genome reference name (hg19, mm10 ...) or data frame from 
   segAnglePo function or data frame from user's mapping data}.
   \item{type}{the type is one of 
      \itemize{
         \item "arc": arcs with variable radii  
         \item "arc2": arcs with the fixed radius
         \item "b": bar charts
         \item "b2": bar charts (opposite side of cutoff value)
         \item "b3": bar charts with the same height
         \item "box": box plots
         \item "chr": plots of chromosomes or segments
         \item "chr2": plots of chromosomes or segments of partial genome
         \item "ci95": 95\% confidence interval lines
         \item "h": histograms
         \item "heatmap": heatmaps
         \item "heatmap2": heatmaps with genomic coordinates
         \item "hightlight.link": link lines for zoom
         \item "hist": polygons for multiple samples
         \item "hl": highlight
         \item "l": lines 
         \item "label": gene labels or text annotations
         \item "label2": gene labels or text annotations with the same circumference coordinate
         \item "lh": horizontal lines
         \item "link.pg": link polygons based on Bezier curve
         \item "link": link lines based on Bezier curve 
         \item "link2": link lines with smaller intra-chromosome arcs
         \item "ls": lines in stair steps
         \item "ml": multiple lines (for more than 1 samples)
         \item "ml2": multiple horizontal lines 
         \item "ml3": multiple lines in stair steps
         \item "ms":  multiple points
         \item "quant75": 75\% quantile lines
         \item "quant90": 90\% quantile lines
         \item "s": dots
         \item "s2": dots with the fixed radius
         \item "s.sd": dots proportional to standard deviation
         \item "ss": dot sizes proportional to the values 
         \item "sv": dot sizes proportional to the variances
       } 
   }
   \item{col.v}{column number. The column value will be displayed. If type=heatmap, the
                number is as the first column.}
   \item{B}{logical: draw background?}
   \item{print.chr.lab}{logical: draw chromosomes or segment labels?}
   \item{col.bar}{logical: draw col.bar? It is for type=heatmap.}
   \item{col.bar.po}{draw col.bar position, e.g. topleft, bottomright.}
   \item{cluster}{logical: cluster and draw Dendrogram at left coner? It is for type=heatmap only.}
   \item{order}{vector: chromosome or segment order}
   \item{scale}{logical: draw scale?}
   \item{cutoff}{numeric: for multiple samples}
   \item{zoom}{vector containing six values: start chromosome, end chromosome, start position, end
              position, start angle and end angle}
   \item{lwd}{numeric: line width}
   \item{cex}{numeric: fond or point sizes}
   \item{col}{character or vector: color names}
   \item{side}{character (in or out): for type=label(2) only}
   \item{...}{...}
}


\author{
 Ying Hu <yhu@mail.nih.gov>
 Chunhua Yan <yanch@mail.nih.gov>
}
\references{
  OmicCircos: an R package for simple and circular visualization of omics data. Cancer Inform. 2014 Jan 16;13:13-20. doi: 10.4137/CIN.S13495. eCollection 2014. PMID: 24526832 [PubMed] PMCID: PMC3921174
}

\examples{

library(OmicCircos);
options(stringsAsFactors = FALSE);

set.seed(1234);

## initial values for simulation data 
seg.num     <- 10;
ind.num     <- 20;
seg.po      <- c(20:50);
link.num    <- 10;
link.pg.num <- 4;
## output simulation data
sim.out <- sim.circos(seg=seg.num, po=seg.po, ind=ind.num, link=link.num, 
  link.pg=link.pg.num);

seg.f     <- sim.out$seg.frame;
seg.v     <- sim.out$seg.mapping;
link.v    <- sim.out$seg.link
link.pg.v <- sim.out$seg.link.pg
seg.num   <- length(unique(seg.f[,1]));

## select segments
seg.name <- paste("chr", 1:seg.num, sep="");
db       <- segAnglePo(seg.f, seg=seg.name);

colors   <- rainbow(seg.num, alpha=0.5);

pdffile  <- "OmicCircos4vignette1.pdf";
pdf(pdffile, 8, 8);
par(mar=c(2, 2, 2, 2));
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");

circos(R=400, cir=db, type="chr",  col=colors, print.chr.lab=TRUE, W=4, scale=TRUE);
circos(R=360, cir=db, W=40, mapping=seg.v, col.v=3, type="l",   B=TRUE, col=colors[1], lwd=2, scale=TRUE);
circos(R=320, cir=db, W=40, mapping=seg.v, col.v=3, type="ls",  B=FALSE, col=colors[9], lwd=2, scale=TRUE);
circos(R=280, cir=db, W=40, mapping=seg.v, col.v=3, type="lh",  B=TRUE, col=colors[7], lwd=2, scale=TRUE);
circos(R=240, cir=db, W=40, mapping=seg.v, col.v=19, type="ml",  B=FALSE, col=colors, lwd=2, scale=TRUE);
circos(R=200, cir=db, W=40, mapping=seg.v, col.v=19, type="ml2", B=TRUE, col=colors, lwd=2);
circos(R=160, cir=db, W=40, mapping=seg.v, col.v=19, type="ml3", B=FALSE, cutoff=5, lwd=2);
circos(R=150, cir=db, W=40, mapping=link.v, type="link", lwd=2, col=colors[c(1,7)]);
circos(R=150, cir=db, W=40, mapping=link.pg.v, type="link.pg", lwd=2, col=sample(colors,link.pg.num));

dev.off()

\dontrun{
demo(OmicCircos4vignette1)
demo(OmicCircos4vignette2)
demo(OmicCircos4vignette3)
demo(OmicCircos4vignette4)
demo(OmicCircos4vignette5)
demo(OmicCircos4vignette6)
demo(OmicCircos4vignette7)
demo(OmicCircos4vignette8)
demo(OmicCircos4vignette9)
demo(OmicCircos4vignette10)
}

}
