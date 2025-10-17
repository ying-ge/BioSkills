\name{fbumLL}
\alias{fbumLL}
\title{
Calculate log likelihood of BUM model
}
\description{
The function calculates the log likelihood of the BUM model.
}
\usage{
fbumLL(parms, x)
}
\arguments{
  \item{parms}{Vector of parameters; lambda and a.
}
  \item{x}{Numerical vector of p-values.
}
}
\value{
Log likelihood.
}
\author{
Marcus Dittrich
}
\examples{
data(pvaluesExample)
pvals <- pvaluesExample[,1]
bum.mle <- fitBumModel(pvals, plot=FALSE)
fbumLL(parms=c(bum.mle$lambda, bum.mle$a), x=pvals)
}