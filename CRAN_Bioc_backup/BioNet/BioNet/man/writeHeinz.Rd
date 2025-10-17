\name{writeHeinz}
\alias{writeHeinz}
\title{
Write input files for HEINZ
}
\description{
Function to write the input files with the node and edge scores for HEINZ. These files are used
to calculate the maximum scoring subnetwork of the graph. The node scores are matched by their names to the nodes of the network, therefore if 
nodes.scores are provided as a vector or matrix, the vector has to be named, respectively the matrix has to be provided
with rownames. If the network contains more nodes than the score vector, the nodes without a score are scored
with the average over all nodes. If the nodes should not be scored and used for the calculation of the 
maximum scoring subnetwork, draw a subnetwork (\code{\link{subNetwork}}) first and use this for the argument network.
The edge scores can be provided as a vector or matrix as the edge.scores argument. If no scores are provided in the arguments, but 
the use.node.scores or use.edge.scores argument is set to TRUE, it will be automatically looked for the "score" attribute of the nodes
and edges of the network.
}
\usage{
writeHeinz(network, file, node.scores=0, edge.scores=0, use.node.score=FALSE, use.edge.score=FALSE)
}
\arguments{
  \item{network}{Network from which to calculate the maximum scoring subnetwork.}
  \item{file}{File to write to.}
  \item{node.scores}{
Numeric vector or matrix of scores for the nodes of the network. Names of the vector or rownames of the matrix have to
correspond to the PPI identifiers of the network. The scores can also be used from the node attribute "score", given one 
score for each node.}
  \item{edge.scores}{Numeric vector of scores for the edges of the network. Edge scores have to be given in the order of the
edges in the network. It is better to append the edge scores as the edge attribute "score" to the network:
\emph{V(network)\$score} and set use.scores to TRUE.}
  \item{use.node.score}{Boolean value, whether to use the node attribute "score" in the network as node scores.}
  \item{use.edge.score}{Boolean value, whether to use the edge attribute "score" in the network as edge scores.}
}
\author{
Daniela Beisser
}
\seealso{
\code{\link{writeHeinzNodes}} and \code{\link{writeHeinzEdges}}
}
\examples{
library(DLBCL)
# use Lymphoma data and graph to find module
data(interactome)
data(dataLym)
# get induced subnetwork for all genes contained on the chip
chipGraph <- subNetwork(dataLym$label, interactome)
score <- dataLym$score001
names(score) <- dataLym$label
\dontrun{writeHeinz(network=chipGraph, file="lymphoma_001", node.scores=score, edge.scores=0)}
}
