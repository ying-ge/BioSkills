\name{scoreOffset}
\alias{scoreOffset}
\title{Change score offset for 2 FDRs}
\description{
Function to change score offset from FDR1 to FDR2.
}
\usage{
scoreOffset(fb, fdr1, fdr2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fb}{
Model from the beta-uniform mixture fitting.
}
  \item{fdr1}{First false discovery rate.}
  \item{fdr2}{Second false discovery rate.}
}
\value{
Offset for the score of the second FDR.
}
\author{Marcus Dittrich}
\seealso{\code{\link{bumOptim}}, \code{\link{fitBumModel}}}
\examples{
data(pvaluesExample)
pvals <- pvaluesExample[,1]
bum <- bumOptim(pvals, starts=10)
scoreOffset(bum, fdr1=0.001, fdr2=0.000001)
}
