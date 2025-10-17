\name{getEdgeList}
\alias{getEdgeList}
\title{Get representation of graph as edgelist}
\description{
A network in \emph{graphNEL} or \emph{igraph} format is converted to an edgelist.
}
\usage{
getEdgeList(network)
}
\arguments{
  \item{network}{Network in \emph{graphNEL} or \emph{igraph} format.}
}
\value{
A matrix whose columns represent the connected edges.}
\author{Marcus Dittrich}
\examples{
library(DLBCL)
data(interactome)
getEdgeList(interactome)[1:10,]
}

