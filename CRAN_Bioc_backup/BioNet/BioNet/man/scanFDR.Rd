\name{scanFDR}
\alias{scanFDR}
\title{Dataframe of scores over a given range of FDRs}
\description{
The function generates a dataframe for a given range of FDRs.
}
\usage{
scanFDR(fb, fdr, labels=names(fb$pvalues))
}
\arguments{
  \item{fb}{Fitted bum model.}
  \item{fdr}{Vector of FDRs.}
  \item{labels}{Data frame labels.}
}
\value{
Dataframe of scores for given p-values and a range of FDRs.
}
\author{Marcus Dittrich}
\seealso{\code{\link{bumOptim}}, \code{\link{fitBumModel}}}
\examples{
data(pvaluesExample)
pvals <- pvaluesExample[,1]
bum <- bumOptim(pvals, starts=10)
scores <- scanFDR(fb=bum, fdr=c(0.1, 0.001, 0.0001))
scores[1:10,]
}