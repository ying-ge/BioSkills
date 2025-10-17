\name{readHeinzGraph}
\Rdversion{1.4}
\alias{readHeinzGraph}
\title{
Convert HEINZ output to graph
}
\description{
Function to convert the HEINZ output to a graph object, or if the output is in matrix form, it will create a list of graphs.
The function needs the node and the original network, from which the module is calculated.
}
\usage{
readHeinzGraph(node.file, network, format=c("graphNEL", "igraph"))
}
\arguments{
  \item{node.file}{
Heinz node output file.
}
  \item{network}{
Original network from which Heinz input was created.
}
  \item{format}{
Graph format of output, either \emph{igraph} or \emph{graphNEL}.
}
}
\value{
Graph object.
}
\author{
Daniela Beisser
}
\examples{
library(DLBCL)
data(interactome)
# precomputed Heinz output files
\dontrun{lib <- file.path(path.package("BioNet"), "extdata") 
module <- readHeinzGraph(node.file=file.path(lib, "lymphoma_nodes_001.txt.0.hnz"), network=interactome, format="graphNEL");
plotModule(module);
}
}

