\name{resamplingPvalues}
\Rdversion{1.4}
\alias{resamplingPvalues}
\title{
Resampling of microarray expression values and test for differential expression.
}
\description{
The function uses a 50\% jackknife resampling to calculate a pseudo-replicate of the expression matrix. The resampled expression values are used thereupon to calculate p-values for the differential expression between the given groups.
Only two-group comparisons are allowed for the perfomed t-test. 
}
\usage{
resamplingPvalues(exprMat, groups, alternative = c("two.sided", "less", "greater"), resampleMat=FALSE)
}
\arguments{
  \item{exprMat}{Matrix with microarray expression values.}
  \item{groups}{Factors for two groups that are tested for differential expression.}
  \item{alternative}{Testing alternatives for the t-test: "two.sided", "less" or "greater".}
  \item{resampleMat}{Boolean value, whether to retrieve the matrix of jacknife resamples or not.}
}
\value{
A result list is returned, consisting of:
\item{p.values}{VNumerical vector of p-values.}
\item{resampleMat}{Matrix of resampled expression values.}
}
\author{
Daniela Beisser
}
\examples{
library(ALL)
data(ALL)
mat <- exprs(ALL)
groups <- factor(c(rep("A", 64), rep("B", 64)))
results <- resamplingPvalues(mat, groups, alternative="greater")
}