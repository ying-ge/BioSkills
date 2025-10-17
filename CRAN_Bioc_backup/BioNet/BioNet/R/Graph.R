
# *********************************************************
# *
# * Graph associated stuff
# *
# *********************************************************



#
# *** faster version of rmSelfLoops 
#
rmSelfLoops <- function(network)
{
  if(is(network, "igraph"))
  {
    return(simplify(network, remove.loops=TRUE))
  }
  if(is(network, "graphNEL"))
  {
    edgelist <- matrix(unlist(strsplit(edgeNames(network), "~")), ncol=numEdges(network), nrow=2)
    pos <- which(edgelist[1,] == edgelist[2,])
    return(removeEdge(from=edgelist[1,pos], to=edgelist[2,pos], network))
  }
}

#
# *** create graph object from source and target vector
#
makeNetwork <- function(source, target, edgemode = "undirected", format = c("graphNEL", "igraph"))
{
  format <- match.arg(format)
  if(format == "graphNEL")
  {
    stopifnot(length(source) == length(target));
    nodes <- unique(c(source,target));
    ge    <- new("graphNEL", nodes = nodes, edgemode=edgemode)
    g     <- addEdge(source, target, ge, 1)
  }
  else
  {
    if(edgemode == "undirected")
    {
      edgemode <- FALSE
    }
    else
    {
      edgemode <- TRUE
    }
    g <- graph.edgelist(cbind(source, target), directed=edgemode)
  }
  return(g)
}

#
# *** Create a subNetwork with matching nodes only
#
.subNetwork0 <- function(nodeList, network)
{
  if(is(network, "igraph"))
  {
     mapping <- seq(1, (length(V(network))))
	 if(is.null(V(network)$name))
	 {
        V(network)$name <- as.character(V(network))
     }
     names(mapping) <- V(network)$name
     nodeList = mapping[nodeList]
     if(any(is.na(nodeList)))
     {
        nodeList = na.omit(nodeList)
        warning("Not all nodes found in network")
     }
     subgr <- induced.subgraph(network, vids=nodeList) 
  }
  else
  {
    subgr <- subGraph(nodes(network)[nodes(network) %in% nodeList], network)
  }
  return(subgr)
}

#
# *** Create a subNetwork and add adjacent nodes from supergraph
#
.subNetwork1 <- function(nodeList, network)
{
  if(is(network, "igraph"))
  {
     mapping <- seq(1, (length(V(network))))
	 if(is.null(V(network)$name))
	 {
        V(network)$name <- as.character(V(network))
     }
     names(mapping) <- V(network)$name
     nodeList = mapping[nodeList]
     if(any(is.na(nodeList)))
     {
        nodeList = na.omit(nodeList)
        warning("Not all nodes found in network")
     }
     neighb <- unique(unlist(neighborhood(network, nodes=nodeList, order=1)))
     subgr <- induced.subgraph(network, vids=neighb) 
  }
  else
  {
    tmp <- unique(c(unlist(adj(network, nodeList)), nodeList))  
    subgr <- subGraph(nodes(network)[nodes(network) %in% tmp], network)
  }
  return(subgr)
}

#
# *** Create a subNetwork 
#
subNetwork <- function(nodeList, network, neighbors=c("none", "first"))
{
  neighbors <- match.arg(neighbors)
  if(neighbors=="first")
  {
    subnet <- .subNetwork1(nodeList, network)
  }
  else{subnet <- .subNetwork0(nodeList, network)}
  return(subnet)
}

#
#  return the largest component as graph
#
largestComp <- function(network)
{
  if(is(network, "graphNEL"))
  {
    cc  <- connectedComp(network)
    idx <- which.max(listLen(cc))
    return(subGraph(cc[[idx]],network))
  }
  else if(is(network, "igraph"))
  {  
      clust <- clusters(network);
      cid   <- which.max(clust$csize);
      lg    <- induced.subgraph(network, V(network)[clust$membership == cid]);
      return(lg); 
    }
}


#
# largest comp with score > level
# -> trivial solution to problem
#
largestScoreComp <- function(network, score, level=0)
{
  if(is(network, "igraph"))
  {
	sorted.score <- na.omit(score[V(network)$name])
	g <- largestComp(induced.subgraph(network, names(sorted.score)[which(sorted.score>level)]))
  }
  if(is(network, "graphNEL"))
  {
	sorted.score <- na.omit(score[nodes(network)])
	g <- largestComp(subGraph(names(sorted.score)[which(sorted.score>level)], network))
  }
  return(g);
}


#
#   returns a network with permuted node labels
#
permutateNodes <- function(network)
{
  if(is(network, "igraph"))
  {
	g <- network;
	V(g)$name <- sample(V(network)$name);
  }
  if(is(network, "graphNEL"))
  {
	g <- network;
	nodes(g) <- sample(nodes(network));
  }
  return(g);
}


# compare.two.networks(network, subnetwork)
# arguments:
#   network
#   subnetwork or 2nd network
#   plot: TRUE or FALSE if cum. degree distribution should be plotted
# values: plot: cumulative degree distribution
#   network parameters:
#     diam.network: diameter of the network
#     diam.subnet: diameter of the subnetwork
#     av.degree.network: average degree of the network
#     av.degree.subnet: average degree of the subnetwork
#     degree.exponent.network: degree exponent gamma of the network
#     degree.exponent.subnet: degree exponent gamma of the subnetwork
#     av.path.length.network: average path length of the network
#     av.path.length.subnet: average path length of the subnetwork
compareNetworks <- function(network1, network2, plot=TRUE)
{
  if(is(network1, "graphNEL"))
  {
    network1 <- igraph.from.graphNEL(network1)
  }
  if(is(network2, "graphNEL"))
  {
    network2 <- igraph.from.graphNEL(network2)
  }
  if(!is.simple(network1) || !is.simple(network2))
  {
    network1 <- largestComp(rmSelfLoops(network1))
    network2 <- largestComp(rmSelfLoops(network2))
    warning("Self-loops were removed and largest component is used for calculation")
  }
  d <- igraph::degree(network1, mode="all")
  d2 <- igraph::degree(network2, mode="all")
  dd <- rev(cumsum(rev(hist(d, -1:max(d), plot=FALSE)$density)))
  dd2 <- rev(cumsum(rev(hist(d2, -1:max(d2), plot=FALSE)$density)))
  degree.exponent.network1 <- power.law.fit(d)
  degree.exponent.network2 <- power.law.fit(d2)
  if(plot == TRUE)
  {
    xmax <- max(length(dd), length(dd2))
    plot(dd, log="xy", xlab="degree", ylab="cumulative frequency", main="Cumulative degree distribution", col=1, xlim=c(1,xmax))
    points(dd2, col=2)
    Lines <- list(bquote("network1, "*gamma==.(round(-degree.exponent.network1$alpha, 2))), bquote("network2, "*gamma==.(round(-degree.exponent.network2$alpha, 2))))
    legend("topright", legend=do.call(expression, Lines), text.col=c(1,2))
  }
  diam.network1 <- diameter(network1, directed=FALSE)
  diam.network2 = diameter(network2, directed=FALSE)
  av.path.length.network1 <- average.path.length(network1, directed=FALSE, unconnected=FALSE)
  av.path.length.network2 <- average.path.length(network2, directed=FALSE, unconnected=FALSE)
  av.degree.network1 <- sum(d)/length(d)
  av.degree.network2 <- sum(d2)/length(d2)
  network.parameters <- list(unname(diam.network1), diam.network2, av.degree.network1, av.degree.network2, -degree.exponent.network1$alpha, -degree.exponent.network2$alpha, av.path.length.network1, av.path.length.network2)
  names(network.parameters) <- c("diam.network1", "diam.network2", "av.degree.network1", "av.degree.network2", "degree.exponent.network1", "degree.exponent.network2", "av.path.length.network1", "av.path.length.network2")
  return(network.parameters)
}

#
#  get graph representation as an edgelist
#
getEdgeList <- function(network)
{
  if(is(network, "graphNEL"))
  {
    em        <- edgeMatrix(network, duplicates=F);
    dat       <- data.frame(cbind(from =I(nodes(network)[em[1,]])));
    dat$from  <- as.character(dat$from);
    dat$to    <- nodes(network)[em[2,]];
    dat$tag1  <- paste(nodes(network)[em[1,]], nodes(network)[em[2,]], sep="~") ;
    dat$tag2  <- paste(nodes(network)[em[2,]], nodes(network)[em[1,]], sep="~") ;
    return(dat);
  }
  if(is(network, "igraph"))
  {
    dat <- get.edgelist(network) 
    return(dat);
  }
}


