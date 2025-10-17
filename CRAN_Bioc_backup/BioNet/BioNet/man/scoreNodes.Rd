\name{scoreNodes}
\alias{scoreNodes}
\title{Score the nodes of a network}
\description{
The function derives scores from the p-values of the nodes of a network.
}
\usage{
scoreNodes(network, fb, fdr=0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{network}{A network in \emph{graphNEL} or \emph{igraph} format.}
  \item{fb}{Fitted bum model.}
  \item{fdr}{False discovery rate.}
}
\value{
Ordered score vector for the nodes of the network.
}
\author{Marcus Dittrich}
\seealso{\code{\link{bumOptim}}, \code{\link{fitBumModel}}}
\examples{
library(DLBCL)
# load p-values
data(dataLym)
# load graph
data(interactome)
# get induced subnetwork for all genes contained on the chip
chipGraph <- subNetwork(dataLym$label, interactome)
p.values <- dataLym$t.pval
names(p.values) <- dataLym$label
bum <- fitBumModel(p.values, plot=TRUE)
scoreNodes(network=chipGraph, fb=bum, fdr=0.001)
}


