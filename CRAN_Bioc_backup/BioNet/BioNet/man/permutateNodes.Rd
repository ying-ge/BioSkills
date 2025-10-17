\name{permutateNodes}
\alias{permutateNodes}
\title{Permute node labels}
\description{
Function to permutate node labels of a given network.
}
\usage{
permutateNodes(network)
}
\arguments{
  \item{network}{Network in \emph{graphNEL} or \emph{igraph} format.}
}
\value{
Network with permutated labels.
}
\author{Marcus Dittrich}
\examples{
library(DLBCL)
data(interactome)
# remove self-loops before permutating the labels
interactome <- rmSelfLoops(interactome)
perm.net <- permutateNodes(interactome)
perm.net
}
