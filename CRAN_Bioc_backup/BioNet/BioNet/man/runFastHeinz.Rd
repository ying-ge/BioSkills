\name{runFastHeinz}
\alias{runFastHeinz}
\title{
Calculate heuristically maximum scoring subnetwork
}
\description{
The function uses an heuristic approach to calculate the maximum scoring subnetwork. Based on the
given network and scores the positive nodes are in the first step aggregated to meta-nodes between which
minimum spanning trees are calculated. In regard to this, shortest paths yield the approximated maximum
scoring subnetwork.
This function can be used if a CPLEX license is not available to calculate the optimal solution. 
}
\usage{
runFastHeinz(network, scores)
}
\arguments{
  \item{network}{
A graph in \emph{igraph} or \emph{graphNEL} format.
}
  \item{scores}{
A named vector, containing the scores for the nodes of the network. All nodes need to be scored in order to run
the algorithm.
}
}
\value{
A subnetwork in the input network format.
}
\author{
Daniela Beisser
}
\seealso{
\code{\link{writeHeinzEdges}}, \code{\link{writeHeinzNodes}}, \code{\link{readHeinzTree}},  \code{\link{readHeinzGraph}}, \code{\link{runHeinz}}
}
\examples{
library(DLBCL)
# load p-values
data(dataLym)
# load graph
data(interactome)
# get induced subnetwork for all genes contained on the chip
interactome <- subNetwork(dataLym$label, interactome)
p.values <- dataLym$t.pval
names(p.values) <- dataLym$label
bum <- fitBumModel(p.values, plot=TRUE)
scores <- scoreNodes(network=interactome, fb=bum, fdr=0.0001)
module <- runFastHeinz(network=interactome, scores=scores)
\dontrun{plotModule(module)}
}
