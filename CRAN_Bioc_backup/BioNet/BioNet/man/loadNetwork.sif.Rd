\name{loadNetwork.sif}
\alias{loadNetwork.sif}
\title{
Load network from Cytoscape sif file
}
\description{
The function loads a network from a Cytoscape sif file. Edge attributes are provided
in the ea.file or vector of ea.files. The node attributes are provided the same way. 
For other formats see \emph{read.graph} in the igraph package.
}
\usage{
loadNetwork.sif(sif.file, na.file, ea.file, format=c("graphNEL", "igraph"), directed=FALSE)
}
\arguments{
  \item{sif.file}{
Cytoscape sif file, containing the network.
}
  \item{na.file}{
File or vector of file with Cytoscape node attibutes. 
}
  \item{ea.file}{
File or vector of file with Cytoscape edge attibutes. 
}
  \item{format}{
Format of output graph, either \emph{graphNEL} or \emph{igraph}.
}
  \item{directed}{
Boolean value for directed or undirected graph.
}
}
\value{
Graph with loaded attributes.
}
\author{
Daniela Beisser
}
\examples{
\dontrun{lib <- file.path(.path.package("BioNet"), "extdata") 
# load interaction file, node attribute file with a node weight of 2 for each node and the edge attribute file with a edge weight of 1 for each edge
network <- loadNetwork.sif(sif.file=file.path(lib,"cytoscape.sif"), na.file=file.path(lib,"n.weight.NA"), ea.file=file.path(lib,"weight.EA"), format="graphNEL", directed=FALSE);
network;
nodeData(network);
edgeData(network);
}
}
