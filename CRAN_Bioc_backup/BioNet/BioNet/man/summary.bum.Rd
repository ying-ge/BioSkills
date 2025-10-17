\name{summary.bum}
\Rdversion{1.4}
\alias{summary.bum}
\title{
Print summary of informations about bum model
}
\description{
The function summarizes information about the bum model.
}
\usage{
\method{summary}{bum}(object, \dots)
}
\arguments{
  \item{object}{
Maximum likelihood estimator object of the beta-uniform mixture fit.}
  \item{\dots}{Other graphic parameters for summary.}
}
\author{
Daniela Beisser
}
\seealso{\code{\link{fitBumModel}}, \code{\link{print.bum}}}
\examples{
data(pvaluesExample)
pvals <- pvaluesExample[,1]
mle <- fitBumModel(pvals, plot=FALSE)
summary(mle)
}