\name{fbum}
\alias{fbum}
\title{Compute the density of the bum distribution}
\description{
Function to compute the density of the beta-uniform mixture model.
}
\usage{
fbum(x, lambda, a)
}
\arguments{
  \item{x}{A numeric value.}
  \item{lambda}{Parameter lambda, mixture parameter, proportion of uniform component}
  \item{a}{Parameter a, shape parameter of beta component}
}
\value{
Value of the density of the bum distribution for x.
}
\references{
S. Pounds, S.W. Morris (2003) Estimating the occurrence of false positives and false negatives in microarray studies by approximating and partitioning the empirical distribution of p-values. \emph{Bioinformatics}, 19(10): 1236-1242.
}
\author{Marcus Dittrich}
\seealso{\code{\link{bumOptim}}, \code{\link{fitBumModel}}}
\examples{
y <- fbum(x=0.5, lambda=0.1, a=0.1)
y
}
