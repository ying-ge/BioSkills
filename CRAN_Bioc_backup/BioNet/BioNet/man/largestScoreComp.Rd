\name{largestScoreComp}
\alias{largestScoreComp}
\title{Component with largest score}
\description{
The function extracts the component of the network with the largest score. All nodes have 
to exceed the given level for the score.
}
\usage{
largestScoreComp(network, score, level=0)
}
\arguments{
  \item{network}{Network in \emph{graphNEL} or \emph{igraph} format.}
  \item{score}{Vector of scores for the network.}
  \item{level}{Cut-off level for the score for the component.}
}
\value{
Subgraph of the network with a score larger than the given level.
}
\author{Marcus Dittrich}
\examples{
library(DLBCL)
data(interactome)
data(dataLym)
network <- rmSelfLoops(interactome)
score <- dataLym$score001
names(score) <- dataLym$label
lComp <- largestScoreComp(network=network, score=score, level=1)
\dontrun{plotModule(lComp)}
}

