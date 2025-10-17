\name{plot.bum}
\Rdversion{1.4}
\alias{plot.bum}
\title{
Quantile-quantile plot for the beta-uniform mixture model
}
\description{
The function plot the theoretical quantiles of the fitted bum model against the
quantiles of the observed p-value distribution.
}
\usage{
\method{plot}{bum}(x, main="QQ-Plot", xlab="Estimated p-value", ylab="Observed p-value", \dots)
}
\arguments{
  \item{x}{
Maximum likelihood estimation object of the fitted bum model.
}
  \item{main}{An overall title for the plot.}
  \item{xlab}{A title for the x axis.}
  \item{ylab}{A title for the y axis.}
  \item{\dots}{Other graphic parameters for the plot.}
}
\seealso{\code{\link{fitBumModel}}, \code{\link{plot.bum}}, \code{\link{bumOptim}}}
\author{
Daniela Beisser
}
\examples{
data(pvaluesExample)
pvals <- pvaluesExample[,1]
mle <- fitBumModel(pvals, plot=FALSE)
plot(mle)
}