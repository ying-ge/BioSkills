\name{pvaluesExample}
\Rdversion{1.4}
\alias{pvaluesExample}
\docType{data}
\title{
Example p-values for aggregation statistics
}
\description{
Data example consisting of a matrix of p-values. Each gene has two corresponding p-values.
These p-values can be aggregated into a p-value of p-values by the method \code{\link{aggrPvals}}.
}
\usage{data(pvaluesExample)}
\examples{
data(pvaluesExample)
pvaluesExample[1:10,]
}

