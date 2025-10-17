\name{plotLLSurface}
\Rdversion{1.4}
\alias{plotLLSurface}
\title{
Log likelihood surface plot
}
\description{
The function plots the log likelihood surface for all a and lambda parameter of the
beta-uniform mixture model.
}
\usage{
plotLLSurface(x, opt=NULL, main="Log-Likelihood Surface", color.palette = heat.colors, nlevels = 32)
}
\arguments{
  \item{x}{Numeric vector of p-values.}
  \item{opt}{List of optimal parameters for a and lambda from the beta-uniform mixture model.}
  \item{main}{The overall title of the plot.}
  \item{color.palette}{Color scheme of the image plot.}
  \item{nlevels}{Number of color levels.}
}
\author{
Marcus Dittrich
}
\examples{
library(DLBCL)
data(dataLym)
pvals <- dataLym$t.pval
names(pvals) <- dataLym$label
mle <- fitBumModel(pvals, plot=FALSE)
plotLLSurface(x=pvals, opt=mle)
}