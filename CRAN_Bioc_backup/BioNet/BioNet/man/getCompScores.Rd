\name{getCompScores}
\alias{getCompScores}
\title{Partition scores for subgraphs of the network}
\description{
The function partitions the scores into scores for each subgraph of the network.
}
\usage{
getCompScores(network, score)
}
\arguments{
  \item{network}{A network in \emph{graphNEL} or \emph{igraph} format.}
  \item{score}{Vector of scores.}
}
\value{
A data frame with the components of the network and the score for each PPI identifier.
}
\author{Marcus Dittrich}
\examples{
library(DLBCL)
data(interactome)
data(dataLym)
# create random subgraph with 100 nodes and their direct neighbors
nodes <- nodes(interactome)[sample(length(nodes(interactome)), 100)]
subnet <- subNetwork(nodeList=nodes, network=interactome, neighbors="first")
score <- dataLym$score001
names(score) <- dataLym$label
getCompScores(score=score, network=subnet)
}
