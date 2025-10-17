\name{makeNetwork}
\alias{makeNetwork}
\title{Create graph from source and target vectors}
\description{
Function to create a graph in \emph{graphNEL} or \emph{igraph} format from a source and 
a target vector.
}
\usage{
makeNetwork(source, target, edgemode="undirected", format=c("graphNEL", "igraph"))
}
\arguments{
  \item{source}{Vector of source nodes.}
  \item{target}{Vector of corresponding target nodes.}
  \item{edgemode}{For an "undirected" or "directed" network.}
  \item{format}{Graph format, eiter \emph{graphNEL} or \emph{igraph}.}
}
\value{
A graph object.
}
\author{Marcus Dittrich}
\seealso{\code{\link{loadNetwork.sif}}, \code{\link{saveNetwork}}}
\examples{
source <- c("a", "b", "c", "d")
target <- c("b", "c", "a", "a")
graph <- makeNetwork(source, target, edgemode="undirected")
}