\name{fdrThreshold}
\alias{fdrThreshold}
\title{
Calculate p-value threshold for given FDR  
}
\description{
The function calculates the p-value threshold tau for a given false discovery rate.
Tau is used for the scoring function.
}
\usage{
fdrThreshold(fdr, fb)
}
\arguments{
  \item{fdr}{False discovery rate.}
  \item{fb}{Model from the beta-uniform mixture fitting.}
}
\value{
P-value threshold tau.
}
\references{ 
S. Pounds, S.W. Morris (2003) Estimating the occurrence of false positives and false negatives in microarray studies by approximating and partitioning the empirical distribution of p-values. \emph{Bioinformatics}, 19(10): 1236-1242.
}
\author{Marcus Dittrich}
\seealso{\code{\link{fbum}}, \code{\link{fitBumModel}} }
\examples{
data(pvaluesExample)
pvals <- pvaluesExample[,1]
bum.mle <- fitBumModel(pvals, plot=FALSE)
tau <- fdrThreshold(fdr=0.001, fb=bum.mle)
tau
}