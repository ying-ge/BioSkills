\name{subNetwork}
\alias{subNetwork}
\title{Create a subGraph}
\description{
The function creates a subgraph with the nodes given in the nodeList or for these nodes including their direct neighbors.
}
\usage{
subNetwork(nodeList, network, neighbors=c("none", "first"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nodeList}{Character vector of nodes, contained in the subgraph.}
  \item{network}{Graph that is used for subgraph extraction.}
  \item{neighbors}{Neighborhood, that is chosen for the subgraph extraction. "none" are only the selected nodes, "first" includes the direct neighbors of the selected nodes.}
}
\value{
A graph object.
}
\author{Marcus Dittrich}
\examples{
library(igraph)
el <- cbind(c("a", "b", "c", "d", "e", "f", "d"), c("b", "c", "d", "e", "f", "a", "b"))
graph <- graph.edgelist(el, directed=TRUE)

node.list <- c("a", "b", "c")
graph2 <- subNetwork(nodeList=node.list, network=graph)
\dontrun{par(mfrow=c(1,2));
plotModule(graph);
plotModule(graph2)}

# or in graphNEL format: 
graph3 <- igraph.to.graphNEL(graph)
graph4 <- subNetwork(nodeList=node.list, network=graph3)
graph3
graph4
}
