\name{readHeinzTree}
\Rdversion{1.4}
\alias{readHeinzTree}
\title{
Convert HEINZ output to tree
}
\description{
Converts the HEINZ output to a tree in graph format. If the output is in matrix form, it will create a list of graphs.
The function needs the node and edge file and the original network from which the module is calculated.
}
\usage{
readHeinzTree(node.file, edge.file, network, format=c("graphNEL", "igraph"))
}
\arguments{
  \item{node.file}{
Heinz node output file.
}
  \item{edge.file}{
Heinz edge output file.
}
  \item{network}{
Original network from which Heinz input was created.
}
  \item{format}{
Output format of the graph, either \emph{igraph} or \emph{graphNEL}.  
}
}
\value{
A graph object.
}
\author{
Daniela Beisser
}
\examples{
library(DLBCL)
data(interactome)
# precomputed Heinz output files
\dontrun{lib <- file.path(.path.package("BioNet"), "extdata") 
module <- readHeinzTree(node.file=file.path(lib, "lymphoma_nodes_001.txt.0.hnz"), edge.file=file.path(lib, "lymphoma_edges_001.txt.0.hnz"), network=interactome, format="graphNEL");
plotModule(module);
}
}
