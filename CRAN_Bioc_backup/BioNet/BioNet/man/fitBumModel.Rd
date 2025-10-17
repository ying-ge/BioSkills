\name{fitBumModel}
\alias{fitBumModel}
\title{
Fit beta-uniform mixture model to a p-value distribution
}
\description{
The function fits a beta-uniform mixture model to a given p-value distribution.
The BUM method was introduced by Stan Pounds and Steve Morris to model the p-value distribution
as a signal-noise decompostion. The signal component is assumed to be B(a,1)-distributed, 
whereas the noise component is uniform-distributed under the null hypothesis.

}
\usage{
fitBumModel(x, plot = TRUE, starts=10)
}
\arguments{
  \item{x}{
Numeric vector of p-values.
}
  \item{plot}{
Boolean value, whether to plot a histogram and qqplot of the p-values with the fitted model.
}
  \item{starts}{
Numeric value giving the number of starts for the optimization.
}
}
\value{
Maximum likelihood estimator object for the fitted bum model.
List of class fb with the following elements:
  \item{lambda}{Fitted parameter \emph{lambda} for the beta-uniform mixture model.}
  \item{a}{Fitted parameter \emph{a} for the beta-uniform mixture model.}
  \item{negLL}{Negative log-likelihood.}
  \item{pvalues}{P-value vector.}
}
\references{
S. Pounds, S.W. Morris (2003) Estimating the occurrence of false positives and false negatives in microarray studies by approximating and partitioning the empirical distribution of p-values. \emph{Bioinformatics}, 19(10): 1236-1242.
}
\author{
Daniela Beisser
}
\examples{
data(pvaluesExample)
pvals <- pvaluesExample[,1]
bum.mle <- fitBumModel(pvals, plot=TRUE)
bum.mle
}