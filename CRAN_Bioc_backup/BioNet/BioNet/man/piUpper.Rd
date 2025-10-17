\name{piUpper}
\alias{piUpper}
\title{Upper bound pi for the fraction of noise}
\description{
The function calculates the upper bound pi for the fraction of noise.
}
\usage{
piUpper(fb)
}
\arguments{
  \item{fb}{Fitted bum model, list with parameters a and lambda.}
}
\value{
Numerical value for the upper bound pi.
}
\author{Marcus Dittrich}
\seealso{\code{\link{bumOptim}}, \code{\link{fitBumModel}}}
\examples{
data(pvaluesExample)
pvals <- pvaluesExample[,1]
bum <- bumOptim(pvals, starts=10)
piUpper(fb=bum)
}