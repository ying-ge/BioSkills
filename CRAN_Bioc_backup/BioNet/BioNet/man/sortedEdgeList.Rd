\name{sortedEdgeList}
\Rdversion{1.4}
\alias{sortedEdgeList}
\title{
Get a sorted edgelist
}
\description{
Function to get a sorted edgelist where the source protein is alphabetically smaller than the target protein from an undirected network.
}
\usage{
sortedEdgeList(network)
}
\arguments{
  \item{network}{
Undirected network in \emph{igraph} or \emph{graphNEL} format.
}
}
\value{
Vector of sorted edges, where the source protein is alphabetically smaller than the target protein.
}
\author{
Daniela Beisser
}
\examples{
library(DLBCL)
data(interactome)
E.list <- sortedEdgeList(interactome)
}

