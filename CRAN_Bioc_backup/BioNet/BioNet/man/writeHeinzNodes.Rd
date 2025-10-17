\name{writeHeinzNodes}
\alias{writeHeinzNodes}
\title{
Write node input file for HEINZ
}
\description{
Function to write an input file with the node scores for HEINZ. This file is used
together with the edge input file to calculate the maximum scoring subnetwork of the 
graph. The scores are matched by their names to the nodes of the network, therefore if 
nodes.scores are provided as a vector or matrix, the vector has to be named, respectively the matrix has to be provided
with rownames.
If the network contains more nodes than the score vector, the nodes without a score are scored
with the average over all nodes. If the nodes should not be scored and used for the calculation of the 
maximum scoring subnetwork, draw a subnetwork \code{\link{subNetwork}} first and use this for the argument network.
}
\usage{
writeHeinzNodes(network, file, node.scores=0, use.score=FALSE)
}
\arguments{
  \item{network}{
Network from which to calculate the maximum scoring subnetwork.
}
  \item{file}{
File to write to.
}
  \item{node.scores}{
Numeric vector or matrix of scores for the nodes of the network. Names of the vector or rownames of the matrix have to
correspond to the PPI identifiers of the network. The scores can also be used from the node attribute "score", given one 
score for each node.
}
  \item{use.score}{
Boolean value, whether to use the node attribute "score" in the network as node scores.
}
}
\details{
Use \code{\link{scoreNodes}} or \code{\link{scoreFunction}}to derive scores from a vector of p-values.
}
\author{
Daniela Beisser
}
\seealso{
 \code{\link{writeHeinzEdges}} and \code{\link{writeHeinz}} 
}
\examples{
#create small network
library(DLBCL)
data(interactome)
small.net <- subNetwork(nodes(interactome)[0:15], interactome)
scores <- c(1:length(nodes(small.net)))
names(scores) <- nodes(small.net)
\dontrun{writeHeinzNodes(network=small.net, file="test_nodes", node.scores=scores)}

# use Lymphoma data and graph to find module
library(DLBCL)
data(interactome)
data(dataLym)
# get induced subnetwork for all genes contained on the chip
chipGraph <- subNetwork(dataLym$label, interactome)
\dontrun{writeHeinzEdges(network=chipGraph, file="lymphoma_edges_001", use.score=FALSE)}
score <- dataLym$score001
names(score) <- dataLym$label
\dontrun{writeHeinzNodes(network=chipGraph, file="lymphoma_nodes_001", node.scores=score)}
}
