\name{scoreFunction}
\alias{scoreFunction}
\title{
 Scoring function for p-values
}
\description{
The function calculates a score for each gene with a given FDR from the fitted beta-uniform mixture model. 
}
\usage{
scoreFunction(fb, fdr=0.01)
}
\arguments{
  \item{fb}{
Model from the beta-uniform mixture fitting.
}
  \item{fdr}{
Numeric constant, from the false discovery rate a p-value threshold is calculated. P-values below this threshold
are considered to be significant and will score positively, p-values a bove the threshold are supposed to arise
from the null model. The FDR can be used to control the size of the maximum scoring subnetwork, by zooming in and out 
in the same region.
}
}
\value{
Score vector for the given p-values.
}
\references{
For details on the score calculation see: M. T. Dittrich, G. W. Klau, A. Rosenwald, T. Dandekar, T. Mueller (2008) Identifying functional modules in protein-protein interaction networks: an integrated exact approach. \emph{(ISMB2008) Bioinformatics}, 24: 13. i223-i231 Jul.
}
\author{
Marcus Dittrich and Daniela Beisser
}
\examples{
data(pvaluesExample)
pvals <- pvaluesExample[,1]
bum.mle <- fitBumModel(pvals, plot=FALSE)
scores <- scoreFunction(fdr=0.1, fb=bum.mle)
scores
}
