\name{plotModule}
\Rdversion{1.4}
\alias{plotModule}
\title{
Plot of the network
}
\description{
The function plots a network from \emph{graphNEL} or \emph{igraph} format, adapted from an igraph plotting function. It is 
just used to visualize the modules. For further plotting options use the plot.igraph function of the igraph package.
The shapes of the nodes can be changed according to the scores argument, then negative scores appear squared.
The color of the nodes can be changed according to the diff.expr argument. Negative values lead to green nodes, positive
values are colored in red. If the vectors are not provided, it will be automatically looked for nodes attributes with the name
\emph{score} and \emph{diff.expr}.
}
\usage{
plotModule(network, layout=layout.fruchterman.reingold, labels=NULL, diff.expr=NULL, scores=NULL, main=NULL, vertex.size=NULL, \dots)
}
\arguments{
  \item{network}{
Network in \emph{graphNEL} or \emph{igraph} format.
}
  \item{layout}{
Layout algorithm, e.g. layout.fruchterman.reingold or layout.kamada.kawai.
}
  \item{labels}{
Labels for the nodes of the network.
}
  \item{diff.expr}{Named numerical vector of differential expression (fold changes) of the nodes in the network. These will be 
  used for coloring of the nodes.  It will be automatically looked for nodes attribute with the name \emph{diff.expr}, if the argument is null.}
  \item{scores}{Named numerical vector of scores of the nodes in the network. These will be used for the shape of the nodes. 
  It will be automatically looked for nodes attribute with the name \emph{score}, if the argument is null.}
  \item{main}{Main title of the plot.}
  \item{vertex.size}{Numerical value or verctor for the size of the vertices.}
  \item{\dots}{Other graphic parameters for the plot.}
}
\author{
Marcus Dittrich and Daniela Beisser
}
\seealso{
 \code{\link{plot3dModule}} 
}
\examples{
library(DLBCL)
data(dataLym)
data(interactome)
interactome <- subNetwork(dataLym$label, interactome)
interactome <- rmSelfLoops(interactome)
fchange <- dataLym$diff
names(fchange) <- dataLym$label
# create random subnetwork
subnet <- largestComp(subNetwork(nodes(interactome)[1:100], interactome))
fchange <- fchange[nodes(subnet)]

# color random subnetwork by the fold change
\dontrun{plotModule(network=subnet, diff.expr=fchange)}
}