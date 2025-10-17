\name{writeHeinzEdges}
\alias{writeHeinzEdges}
\title{
Write edge input file for HEINZ
}
\description{
Function to write an input file for HEINZ with edge scores. If no edge scores are used, they are set to 0.
In order to run HEINZ, a node input and edge input file are needed.
}
\usage{
writeHeinzEdges(network, file, edge.scores=0, use.score=FALSE)
}
\arguments{
  \item{network}{
Network from which to calculate the maximum scoring subnetwork.
}
  \item{file}{
File to write to.
}
  \item{edge.scores}{
Numeric vector of scores for the edges of the network. Edge scores have to be given in the order of the
edges in the network. It is better to append the edge scores as the edge attribute "score" to the network:
\emph{V(network)\$score} and set use.score to TRUE.
}
  \item{use.score}{
Boolean value, whether to use the edge attribute "score" in the network as edge scores.
}
}
\seealso{
\code{\link{writeHeinzNodes}} and \code{\link{writeHeinz}}
}
\author{Daniela Beisser}
\examples{
library(DLBCL)
# use Lymphoma data and graph to find module
data(interactome)
data(dataLym)
# get induced subnetwork for all genes contained on the chip
chipGraph <- subNetwork(dataLym$label, interactome)
# remove self loops
graph <- rmSelfLoops(chipGraph)
\dontrun{writeHeinzEdges(network=graph, file="lymphoma_edges_001", use.score=FALSE)}
score <- dataLym$score001
names(score) <- dataLym$label
\dontrun{writeHeinzNodes(network=graph, file="lymphoma_nodes_001", node.scores=score)}


# write another edge file with edge scores 
library(igraph)
data(interactome)
interactome <- igraph.from.graphNEL(interactome)
small.net <- subNetwork(V(interactome)[1:16]$name, interactome)
scores <- c(1:length(E(small.net)))
E(small.net)$score <- scores
\dontrun{writeHeinzEdges(network=small.net, file="test_edges", use.score=TRUE)}
}
