\name{plot3dModule}
\Rdversion{1.4}
\alias{plot3dModule}
\title{
3D plot of the network
}
\description{
The function plots a network from \emph{graphNEL} or \emph{igraph} format in 3D using a modified function from the package igraph and requires the package
rgl which uses openGL. The 3D plot can be zoomed, rotated, shifted on the canvas.  This function is just used to visualize the modules. For further plotting options
use the rglplot function of the igraph package. If a score attribute is provided in the graph this will be used for the
coloring of the nodes. Otherwise a vector of values can be given by the diff.or.score argument. The vector has to contain positive and negative values,
either scores or values for differential expression (fold changes). Labels for the nodes can be provided by the labels argument, otherwise it will be automatically
looked for a \emph{geneSymbol} attribute of the nodes.
}
\usage{
plot3dModule(network, labels=NULL, windowSize = c(100,100,1500,1000), diff.or.scores=NULL, red=c("negative", "positive"), \dots)
}
\arguments{
  \item{network}{
Network in \emph{graphNEL} or \emph{igraph} format.
}
  \item{labels}{
Labels for the nodes of the network. Otherwise it will be automatically looked for a \emph{geneSymbol} attribute of the nodes.
}
  \item{windowSize}{Numerical vector of size four to set the size of the rgl device.}
  \item{diff.or.scores}{Named numerical vector of differential expression (fold changes) or scores of the nodes in the network. These will be used for node coloring. Otherwise a \emph{score} attribute of the nodes will be automatically used.}
  \item{red}{Either "negative" or "positive", to specify which values are to be colored red in the plot.}
  \item{\dots}{Other graphic parameters for the plot.}
}
\author{
Daniela Beisser
}
\seealso{\code{\link{save3dModule}}, \code{\link{plotModule}}}
\examples{
library(DLBCL)
data(interactome)
data(dataLym)
interactome <- subNetwork(dataLym$label, interactome)
interactome <- rmSelfLoops(interactome)
fchange <- dataLym$diff
names(fchange) <- dataLym$label
subnet <- largestComp(subNetwork(nodes(interactome)[1:100], interactome))
diff <- fchange[nodes(subnet)]

\dontrun{library(rgl);
plot3dModule(network=subnet, diff.or.score=diff)}
}