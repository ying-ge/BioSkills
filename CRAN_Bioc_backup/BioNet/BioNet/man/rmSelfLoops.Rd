\name{rmSelfLoops}
\alias{rmSelfLoops}
\title{Remove self-loops in a graph}
\description{
The function removes self-loops, edges that start and end in the same node, from the network.
}
\usage{
rmSelfLoops(network)
}
\arguments{
  \item{network}{A graph object, either in \emph{graphNEL} or \emph{igraph} format.}
}
\value{
The graph with the removed edges.
}
\author{Marcus Dittrich}
\examples{
graph <- makeNetwork(c("a","b","c","d","e","a"), c("b","c","d","e","e","e"))
graph2 <- rmSelfLoops(graph)
edges(graph)
edges(graph2)
}