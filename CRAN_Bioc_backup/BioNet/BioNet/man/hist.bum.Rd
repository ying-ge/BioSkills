\name{hist.bum}
\Rdversion{1.4}
\alias{hist.bum}
\title{
Histogram of the p-value distribution with the fitted bum model
}
\description{
The function plots a histogram of the p-values together with the fitted bum-model.
}
\usage{
\method{hist}{bum}(x, breaks=50, main="Histogram of p-values", xlab="P-values", ylab="Density", \dots)
}
\arguments{
  \item{x}{
Maximum likelihood estimator object of the beta-uniform mixture fit.}
  \item{breaks}{Breaks for the histogram.}
  \item{main}{An overall title for the plot.}
  \item{xlab}{A title for the x axis.}
  \item{ylab}{A title for the y axis.}
  \item{\dots}{Other graphic parameters for the plot.}
}
\author{
Daniela Beisser
}
\seealso{\code{\link{fitBumModel}}, \code{\link{hist.bum}}, \code{\link{bumOptim}}}
\examples{
data(pvaluesExample)
pvals <- pvaluesExample[,1]
mle <- fitBumModel(pvals, plot=FALSE)
hist(mle)
}
