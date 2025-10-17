\name{saveNetwork}
\alias{saveNetwork}
\title{
Save undirected network in various formats
}
\description{
The function saves a graph in a Cytoscape readable format: either in XGMML format, or
as two tables, one for the nodes with attributes and one for the edges with attributes, or as .sif file.
Or other standard formats like tab separated, .tgf, .net
}
\usage{
saveNetwork(network, name="network", file, type=c("table", "XGMML", "sif", "tab", "tgf", "net"))
}
\arguments{
  \item{network}{
Network to save.
}
  \item{name}{
Name of the network, only needed for the XGMML format.
}
  \item{file}{
File to save to.
}
  \item{type}{
Type in which graph shall be saved.
}
}
\details{
The format types are "XGMML", "table", "sif", "tab", "tgf" and "net".
XGMML (eXtensible Graph Markup and Modeling Language) is an 
XML format based on GML which is used for graph description. Edges, nodes and their 
affiliated attributes are all saved in one file. 
In the table format two tables are created, one for the nodes with attributes and one for the edges with attributes.
The sif format creates a .sif file for the network and an node attribute (.NA) or edge attribute (.EA) for 
each attribute. The name of the attribute is the filename.  
Tab writes only the edges of the network in a tabular format. Tgf save the network to simple .tgf format.
The net format writes a Pajek readable file of the network and the ET type saves the edge tags to file.
}
\author{
Daniela Beisser and Marcus Dittrich
}
\examples{
library(DLBCL)
#create small network
library(igraph)
data(interactome) 
interactome <- igraph.from.graphNEL(interactome)
small.net <- subNetwork(V(interactome)[1:16]$name, interactome)
E(small.net)$e.weight <- rep(1,length(E(small.net)))
V(small.net)$n.weight <- rep(2,length(V(small.net)))
summary(small.net)
\dontrun{saveNetwork(small.net, file="example_network", name="small.net", type="XGMML")}
}

