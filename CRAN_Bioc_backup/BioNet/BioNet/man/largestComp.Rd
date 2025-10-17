\name{largestComp}
\alias{largestComp}
\title{Extract largest component of network}
\description{
The function extracts the largest component of a network.
}
\usage{
largestComp(network)
}
\arguments{
  \item{network}{A graph in \emph{graphNEL} or \emph{igraph} format.}
}
\value{
A new graph object that represents the largest component of the given network.
}
\author{Marcus Dittrich}
\examples{
library(DLBCL)
data(interactome)
interactome
largestComp(interactome)
}
