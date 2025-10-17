\name{bumOptim}
\alias{bumOptim}
\title{ 
Fitting a beta-uniform mixture model to p-value distribution
}
\description{
The function fits a beta-uniform mixture model to a given p-value distribution.
}
\usage{
bumOptim(x, starts=1, labels=NULL)
}
\arguments{
  \item{x}{
Numerical vector of p-values, has to be named with the gene names or the gene names can be given in the labels paramater.  
  }
  \item{starts}{
Number of start points for the optimization.
}
  \item{labels}{
Gene names for the p-values.
}
}
\value{
List of class fb with the following elements:
  \item{lambda}{Fitted parameter \emph{lambda} for the beta-uniform mixture model.}
  \item{a}{Fitted parameter \emph{a} for the beta-uniform mixture model.}
  \item{negLL}{Negative log-likelihood.}
  \item{pvalues}{P-value vector.}
}
\references{
M. T. Dittrich, G. W. Klau, A. Rosenwald, T. Dandekar, T. Mueller (2008) Identifying functional modules in protein-protein interaction networks: an integrated exact approach. \emph{(ISMB2008) Bioinformatics}, 24: 13. i223-i231 Jul.

S. Pounds, S.W. Morris (2003) Estimating the occurrence of false positives and false negatives in microarray studies by approximating and partitioning the empirical distribution of p-values. \emph{Bioinformatics}, 19(10): 1236-1242.
}
\author{Marcus Dittrich and Daniela Beisser}
\seealso{\code{\link{fitBumModel}}, \code{\link{plot.bum}}, \code{\link{hist.bum}}}
\examples{
data(pvaluesExample)
pvals <- pvaluesExample[,1]
bum <- bumOptim(x=pvals, starts=10)
bum
}
