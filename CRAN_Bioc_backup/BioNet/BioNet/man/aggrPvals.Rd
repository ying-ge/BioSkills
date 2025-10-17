\name{aggrPvals}
\alias{aggrPvals}
\title{
Aggregate several p-values into one p-value
}
\description{
The function aggregates several p-values into one p-value of p-values based on the order 
statistics of p-values. An overall p-value is given by the \emph{i}th order statistic.
}
\usage{
aggrPvals(pval.matrix, order, plot=TRUE)
}
\arguments{
  \item{pval.matrix}{
Numeric matrix of p-values, columns represent different sets of p-values
}
  \item{order}{
Numeric constant, the order statistic that is used for the aggregation.
}
  \item{plot}{
Boolean value whether to plot p-value distributions.
}
}
\value{
Aggregated p-value of the given order.
}
\author{
Daniela Beisser
}
\examples{
data(pvaluesExample)
aggrPvals(pval.matrix=pvaluesExample, order=2)
}
