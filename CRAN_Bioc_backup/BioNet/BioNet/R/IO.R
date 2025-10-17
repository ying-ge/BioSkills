
# *********************************************************
# *
# * IO-Routines
# *
# *********************************************************



# *********************************************************
# * Graph IO
# *********************************************************

# save.network(network, name="network", file, type="table")
# arguments:
#   network: network to save
#   name: name of the network for cytoscape
#   file: fiel to save to
#   type: table or XGMML format
saveNetwork <- function(network, name="network", file, type=c("table", "XGMML", "sif", "tab", "tgf", "net"))
{
  file <- .cleanFile(file)
  type <- match.arg(type)
  if(type == "XGMML")
  {
    requireNamespace("XML")
    addNode <- XML::addNode
    if(is(network, "graphNEL"))
    {
      network <- igraph.from.graphNEL(network)  
    }
    top <- .XGMML.destription(name=name)
    print("...adding nodes")
    # append nodes
    nodes <- .XGMML.nodes(network=network)
    top <- XML::append.xmlNode(top, nodes) 
    print("...adding edges")
    # append edges
    edges <- .XGMML.edges(network=network)
    top <- XML::append.xmlNode(top, edges) 
    # save to file as xgmml
    print("...writing to file")
    XML::saveXML(top, file=paste(file, ".xgmml", sep=""), encoding="UTF-8")
	if("package:XML" %in% search()){detach("package:XML")}
	addNode <- graph::addNode
  }
  if(type == "table")
  {
    if(is(network, "graphNEL"))
    {
      network <- igraph.from.graphNEL(network)  
    }
    .graph.table(network=network, file=file)
  }
  if(type == "sif")
  {
    if(is(network, "graphNEL"))
    {
      network <- igraph.from.graphNEL(network)  
    }
    edges <- .graph.sif(network=network, file=file)
	if(length(list.edge.attributes(network))!=0)
	{
		.graph.eda(network=network, file=file, edgelist.names=edges)
	}
	if(length(list.vertex.attributes(network))!=0)
	{
		.graph.noa(network=network, file=file)
	}
  }
  if(type == "tab")
  {
    if(is(network, "igraph"))
    {
		nE <- ecount(network)
		network <- simplify(network, remove.multiple = TRUE)
		if(nE!=ecount(network))
		{
			warning("Multiple edges are not allowed for the graphNEL format, they had to be removed")
		}
		network <- igraph.to.graphNEL(network)
		.saveGraph.tab(graph=network, filename=file)
    }  
  }
  if(type == "tgf")
  {
		if(is(network, "igraph"))
		{
		nE <- ecount(network)
		network <- simplify(network, remove.multiple = TRUE)
		if(nE!=ecount(network))
		{
			warning("Multiple edges are not allowed for the graphNEL format, they had to be removed")
		}
		network <- igraph.to.graphNEL(network)
		.saveGraph.tgf(graph=network, filename=file) 
    }
  }
  if(type == "net")
  {
    if(is(network, "igraph"))
    {
		nE <- ecount(network)
		network <- simplify(network, remove.multiple = TRUE)
		if(nE!=ecount(network))
		{
			warning("Multiple edges are not allowed for the graphNEL format, they had to be removed")
		}
		network <- igraph.to.graphNEL(network)
		.saveGraph.net(graph=network, filename=file) 
    }  
  }
}

# load.network(sif.file, na.file, ea.file)
# arguments:
#   sif.file: cytoscape sip file
#   na.file: cytoscape node attribute file (vector)
#   ea.file: cytoscape edge attribute file (vector)
#   format: "igraph" or "graphNEL"
#   directed: TRUE or FALSE
# values: igraph network
loadNetwork.sif <- function(sif.file, na.file=NULL, ea.file=NULL, format=c("graphNEL", "igraph"), directed=FALSE)
{
  format <- match.arg(format)
  print("...loading network")
  network <- read.table(file=sif.file, colClasses = "character")
  network <- graph.edgelist(as.matrix(network[,-2]), directed=directed)
  if(!is.null(na.file))
  {
    print("...loading node attributes")
    network <- .add.node.attrs(network, na.file)
  }
  if(!is.null(ea.file))
  {
    print("...loading edge attributes")
    network <- .add.edge.attrs(network, ea.file)
  }
  if(format == "graphNEL")
  {
	nE <- ecount(network)
	network <- simplify(network, remove.multiple = TRUE)
	if(nE!=ecount(network))
	{
		warning("Multiple edges are not allowed for the graphNEL format, they had to be removed")
	}
    network <- igraph.to.graphNEL(network)
  }
  return(network)
} 

#
# *** read a Graph from tab delimited edgelist
#
loadNetwork.tab <- function(file, header=TRUE, directed=FALSE, format=c("graphNEL", "igraph"))
{
  format <- match.arg(format)
  if(directed)
  {
    directed <- "directed"
  }
  if(!directed)
  {
    directed <- "undirected"
  }
  if(header) sk=1 else sk=0
  links <- scan(file, what=list('character', 'character'), sep="\t", skip=sk)
  nodes <- unique(c(links[[1]], links[[2]]))
  ge    <- new("graphNEL", nodes=nodes, edgemode=directed)
  graph     <- addEdge(links[[1]], links[[2]], ge, 1)
  if(format == "igraph")
  {
    graph <- igraph.from.graphNEL(graph)
  }
  return(graph)
}

# internal methods #############################################################

# internal method to add node attributes
.add.node.attrs <- function(network, na.file)
{
  for(i in 1:length(na.file))
  {
    # first line
    attr <- strsplit(readLines(con=na.file[i], n=1)," ")[[1]][1]
    na.input <- read.table(file=na.file[i], skip=1, sep="=", colClasses = "character")
    node.attr <- sapply(as.character(na.input[,2]), function(t){substr(t, 2, nchar(t))})   
    names(node.attr) <- sapply(as.character(na.input[,1]), function(t){substr(t, 1, nchar(t)-1)})
    network <- set.vertex.attribute(network, index=V(network), name=attr, value=as.vector(node.attr[V(network)$name]))    
  }
  return(network) 
}

# internal function to add edge attributes
.add.edge.attrs <- function(network, ea.file)
{
  for(i in 1:length(ea.file))
  {
    # first line
    attr <- strsplit(readLines(con=ea.file[i], n=1)," ")[[1]][1]
    ea.input <- read.table(file=ea.file[i], skip=1, sep="=", colClasses = "character")
    interac.id <- sapply(as.character(ea.input[,1]), function(t){substr(t, 1, nchar(t)-1)})
    interac.id <- sub("\\s\\([a-z\\s0-9]+\\)\\s", " ", interac.id)
    edge.attr <- sapply(as.character(ea.input[,2]), function(t){substr(t, 2, nchar(t))})
    names(edge.attr) <- interac.id  
    edgelist <- get.edgelist(network)   
    edgelist1 <- paste(edgelist[,1], edgelist[,2], sep=" ")  
    edgelist2 <- paste(edgelist[,2], edgelist[,1], sep=" ")  
    value1 <- as.vector(edge.attr[edgelist1])
    value2 <- as.vector(edge.attr[edgelist2])
    value1[which(is.na(value1))] <- value2[which(is.na(value1))]
    network <- set.edge.attribute(network, index=E(network), name=attr, value=value1)  
  }
  return(network)
}

# internal method to save the graph in table format
.graph.table <- function(network, file)
{
  # create vertex attributes
  attrib <- list.vertex.attributes(network)
  if(length(attrib)!=0)
  {
    node.attribs <- matrix(data=NA, ncol=length(attrib), nrow=length(V(network)))
    for(i in 1:length(attrib))
    {
      node.attribs[,i] <- get.vertex.attribute(network, attrib[i])
    }
    node.attribs <- cbind(V(network)$name, node.attribs) 
    colnames(node.attribs) <- c("id", attrib)
    write.table(node.attribs, file=paste(file, "_n.txt", sep=""), row.names=FALSE, sep="\t")
  }
  attrib <- list.edge.attributes(network)
  if(length(attrib)!=0)
  {
    edge.attribs <- matrix(data=NA, ncol=length(attrib), nrow=length(E(network)))
    for(i in 1:length(attrib))
    {
      edge.attribs[,i] <- get.edge.attribute(network, attrib[i])
    }
    edgelist.names <- get.edgelist(network, names=TRUE)
    edge.attribs <- cbind(edgelist.names, edge.attribs)
    colnames(edge.attribs) <- c("nodeA", "nodeB", attrib)
    write.table(edge.attribs, file=paste(file, "_e.txt", sep=""), row.names=FALSE, sep="\t")
  }
}    

# internal function to write cytoscape .sif file
.graph.sif <- function(network, file)
{
  edgelist.names <- get.edgelist(network, names=TRUE)
  edgelist.names <- cbind(edgelist.names[,1], rep("pp", length(E(network))), edgelist.names[,2])
  write.table(edgelist.names, row.names=FALSE, col.names=FALSE, file=paste(file, ".sif", sep=""), sep="\t", quote=FALSE)
  return(edgelist.names) 
} 

# internal method to write cytoscape node attribute files
.graph.noa <- function(network, file)
{
  attrib <- list.vertex.attributes(network)
  for(i in 1:length(attrib))
  {
    if(is(get.vertex.attribute(network, attrib[i]))[1] == "character")
    {
      type <- "String"
    }
    if(is(get.vertex.attribute(network, attrib[i]))[1] == "integer")
    {
      type <- "Integer"
    }
    if(is(get.vertex.attribute(network, attrib[i]))[1] == "numeric")
    {
      type <- "Double"
    }
    noa <- cbind(V(network)$name, rep("=", length(V(network))), get.vertex.attribute(network, attrib[i]))
    first.line <- paste(attrib[i], " (class=java.lang.", type, ")", sep="")
    write(first.line, file=paste(file, "_", attrib[i], ".NA", sep=""), ncolumns = 1, append=FALSE, sep=" ")
    write.table(noa, row.names = FALSE, col.names = FALSE, file=paste(file, "_", attrib[i], ".NA", sep=""), sep=" ", append=TRUE, quote=FALSE)
  }
}

# internal method to write cytoscape edge attribute files
.graph.eda <- function(network, file, edgelist.names)
{
  attrib <- list.edge.attributes(network)
  for(i in 1:length(attrib))
  {
    if(is(get.edge.attribute(network, attrib[i]))[1] == "character")
    {
      type <- "String"
    }
    if(is(get.edge.attribute(network, attrib[i]))[1] == "integer")
    {
      type <- "Integer"
    }
    if(is(get.edge.attribute(network, attrib[i]))[1] == "numeric")
    {
      type <- "Double"
    }
    eda <- cbind(cbind(edgelist.names[,1], rep("(pp)", length(E(network))), edgelist.names[,3]), rep("=", length(E(network))), get.edge.attribute(network, attrib[i]))
    first.line <- paste(attrib[i], " (class=java.lang.", type, ")", sep="")
    write(first.line, file=paste(file, "_", attrib[i], ".EA", sep=""), ncolumns=1, append=FALSE, sep =" ")
    write.table(eda, row.names = FALSE, col.names = FALSE, file=paste(file, "_", attrib[i], ".EA", sep=""), sep=" ", append=TRUE, quote=FALSE)
  }
}
                                                    
# internal method to create the first part of the XGMML-file, description
.XGMML.destription <- function(name)
{
  requireNamespace("XML")
  # create top node
  # top part of xml
  top <- XML::xmlNode("graph", attrs = c(label=name, "xmlns:dc"="http://purl.org/dc/elements/1.1/", "xmlns:xlink"="http://www.w3.org/1999/xlink", "xmlns:rdf"="http://www.w3.org/1999/02/22-rdf-syntax-ns#", "xmlns:cy"="http://www.cytoscape.org", xmlns="http://www.cs.rpi.edu/XGMML"))
  top <- XML::append.xmlNode(top, XML::xmlNode("att", attrs=c(name="documentVersion", value="1.1")))
  
  d <- XML::xmlNode("rdf:Description", attrs=c("rdf:about"="http://www.cytoscape.org/"))
  d <- XML::append.xmlNode(d,  XML::xmlNode("dc:type", "Protein-Protein Interaction"))
  d <- XML::append.xmlNode(d,  XML::xmlNode("dc:description", "N/A"))
  d <- XML::append.xmlNode(d,  XML::xmlNode("dc:identifier", "N/A"))
  d <- XML::append.xmlNode(d,  XML::xmlNode("dc:date", Sys.time()))
  d <- XML::append.xmlNode(d,  XML::xmlNode("dc:title", name))
  d <- XML::append.xmlNode(d,  XML::xmlNode("dc:format", "BioNet-Cytoscape-XGMML"))
  
  c <- XML::xmlNode("att", attrs=c(name="networkMetadata"), XML::xmlNode("rdf:RDF", d))
  top <- XML::append.xmlNode(top, c)
  return(top)
}

# internal method for the addition of nodes to xml
.XGMML.nodes <- function(network)
{
  requireNamespace("XML")
  # create node-nodes
  c.node <- rep("node", length(V(network)))
  nodes <- lapply(c.node, XML::xmlNode)
  
  # create node attributes
  attrib <- list.vertex.attributes(network)
  node.attribs <- matrix(data=NA, nrow=length(attrib), ncol=length(V(network)))
  for(i in 1:length(attrib))
  {
      if(is(get.vertex.attribute(network, attrib[i]))[1] == "character")
      {
        type <- "string"
      }
      if(is(get.vertex.attribute(network, attrib[i]))[1] == "integer")
      {
        type <- "integer"
      }
      if(is(get.vertex.attribute(network, attrib[i]))[1] == "numeric")
      {
        type <- "real"
      }
      node.attribs[i,] =paste("att type=", "\"", type, "\"", " name=", "\"", attrib[i], "\"", " value=", "\"", get.vertex.attribute(network, attrib[i]), "\"", sep="")
  }
  node.attribs <- matrix(lapply(node.attribs, XML::xmlNode), nrow = length(attrib), ncol = length(V(network)))
  if(is.null(V(network)$name))
  {
    V(network)$name <- as.character(V(network))
  }
  node.label <- V(network)$name
  node.id <- as.vector(V(network))
  
  # append node attributes
  for(i in 1:length(V(network)))
  {
    nodes[[i]] <- XML::addAttributes(nodes[[i]], label = node.label[i], id=node.id[i])
    nodes[[i]] <- XML::append.xmlNode(nodes[[i]], node.attribs[,i])
  }
  
  return(nodes)
}

# internal method for the addition of edges to XGMML
.XGMML.edges <- function(network)
{
  requireNamespace("XML")
  # create edge-nodes
  c.edge <- rep("edge", length(E(network)))
  edges <- lapply(c.edge, XML::xmlNode)
  
  edgelist.names <- get.edgelist(network, names=TRUE)
  edgelist.names <- paste(edgelist.names[,1], edgelist.names[,2], sep=" (pp) ")
  edgelist.ids <- get.edgelist(network, names=FALSE)
   
  # create edge attributes
  attrib <- list.edge.attributes(network)
  edge.attribs <- matrix(data=NA, nrow=length(attrib), ncol=length(E(network)))
  for(i in 1:length(attrib))
  {
      if(is(get.edge.attribute(network, attrib[i]))[1] == "character")
      {
        type <- "string"
      }
      if(is(get.edge.attribute(network, attrib[i]))[1] == "integer")
      {
        type <- "integer"
      }
      if(is(get.edge.attribute(network, attrib[i]))[1] == "numeric")
      {
        type <- "real"
      }
      edge.attribs[i,]  <- paste("att type=", "\"", type, "\"", " name=", "\"", attrib[i], "\"", " value=", "\"", get.edge.attribute(network, attrib[i]), "\"", sep="")
  }
  edge.attribs <- matrix(lapply(edge.attribs, XML::xmlNode), nrow=length(attrib), ncol=length(E(network)))
  
  # append edge attributes
  for(i in 1:length(E(network)))
  {
    edges[[i]] <- XML::addAttributes(edges[[i]], label=edgelist.names[i], source=edgelist.ids[i,1], target=edgelist.ids[i,2])
    edges[[i]] <- XML::append.xmlNode(edges[[i]], c(edge.attribs[,i]))
  }
  
  return(edges)
}

#
# *** save graph to simple .tgf format
#
.saveGraph.tgf <- function(graph, filename)
{
  write.table(nodes(graph), file=filename, col.names=F, quote=F)
  write("\n#\n", append=T, file=filename)
  ew <- eWV(graph, edgeMatrix(graph), sep=" ")
  write(names(ew), file=paste(filename, ".tgf", sep=""), append=T)
}

#
# *** save Graph in .net (Pajek) format
#
.saveGraph.net <- function(graph, filename)
{
  n <- length(nodes(graph));
  write(paste("*Vertices ", n), file=filename);
  write.table(dQuote(nodes(graph)), file=filename, col.names=F, quote=F, append=T)
  if(isDirected(graph))
  {
    write("*Arcs", append=T, file=filename);
  }
  else write("*Edges", append=T, file=filename);

  ew <- eWV(graph, edgeMatrix(graph), sep=" ")
  write(names(ew), file=paste(filename, ".net", sep=""), append=T)
}

#
# *** 
#
.saveGraph.tab <- function(graph, filename)
{
  en <- edgeNames(graph);
  x  <- gsub("\\~", "\t", en);
  cat(x, file=paste(filename, ".txt", sep=""), sep="\n");
}


.cleanFile = function(file)
{
  file.vector <- unlist(strsplit(file, "\\."))
  if(file.vector[length(file.vector)] %in% c("txt", "sif", "tab", "XGMML", "tgf", "NOA", "EDA", "xgmml", "eda", "noa", "net", "pdf"))
  {
    file.vector <- file.vector[-length(file.vector)]
    file <- paste(file.vector, collapse=".")
  }
  return(file) 
}


################################################################################
