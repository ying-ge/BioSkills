\name{print.bum}
\Rdversion{1.4}
\alias{print.bum}
\title{
Print information about bum model
}
\description{
The function prints information about the bum model.
}
\usage{
\method{print}{bum}(x, \dots)
}
\arguments{
  \item{x}{
Maximum likelihood estimator object of the beta-uniform mixture fit.}
  \item{\dots}{Other graphic parameters for print.}
}
\author{
Marcus Dittrich
}
\seealso{\code{\link{fitBumModel}}, \code{\link{summary.bum}}}
\examples{
data(pvaluesExample)
pvals <- pvaluesExample[,1]
mle <- fitBumModel(pvals, plot=FALSE)
print(mle)
}