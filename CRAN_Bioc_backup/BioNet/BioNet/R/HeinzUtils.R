
# *********************************************************
# *
# * Heinz utilities for version > 1.6
# *
# *********************************************************

# write file that can be read from heinz as edge input file
# write.heinz.e.input(network, file.out, edge.scores=rep(0, length(E(network))), use.weights = TRUE)
# arguments:
#   network
#   file: a character string naming the file to write to
#   edge.scores: vector of edge scores belonging to the network, else zeros
#   use.weights: boolean, if the edges of the network are weighted, the edge weights can be used as edge scores
writeHeinzEdges <- function(network, file, edge.scores=0, use.score=FALSE)
{
  file <- .cleanFile(file)
  if(is(network, "graphNEL"))
  {
    network <- igraph.from.graphNEL(network)
  }
  if(length(edge.scores) == 1 && edge.scores == 0)
  {
    edge.scores <- rep(0, length(E(network)))
  }
  if(use.score)
  {
    if(!is.null(E(network)$score))
    {
      edge.scores <- E(network)$score
    }
  }
  edges <- get.edgelist(network)
  table.edges <- cbind(edges, edge.scores)
  c.name=c("#nodeA", "nodeB")
  for(i in 3:dim(table.edges)[2])
  {
    expect <- sum(as.vector(table.edges[, i], mode="numeric"), na.rm = TRUE)/length(table.edges[,i])
    table.edges[which(is.na(table.edges[,i])),i] <- expect
    c.name <- c(c.name, paste("score", i-2, sep=""))
  }
  colnames(table.edges) <- c.name
  write.table(table.edges, file=paste(file, ".txt", sep=""), append=FALSE, quote=FALSE, sep="\t", na="NA", dec=".", row.names=FALSE, col.names=TRUE)
  return(TRUE)
}

# write file that can be read from heinz as node input file
# write.heinz.n.input(network, file.out, edge.scores=rep(0, length(E(network))), use.weights = TRUE)
# arguments:
#   network
#   file: a character string naming the file to write to
#   node.scores: vector of node scores belonging to the network, else zeros
writeHeinzNodes <- function(network, file, node.scores=0, use.score=FALSE)
{
  file <- .cleanFile(file)
  if(is(network, "graphNEL"))
  {
    network <- igraph.from.graphNEL(network)
  }
  if(length(node.scores) == 1 && node.scores == 0)
  {
    node.scores <- rep(0, length(V(network)))
  }
  if(is.null(V(network)$name))
  {
    V(network)$name <- as.character(V(network))
  }
  if(use.score)
  {
    if(!is.null(V(network)$score))
    {
      node.scores <- V(network)$score
      names(node.scores) <- V(network)$name
    }
  }
  if(is(node.scores, "data.frame"))
  {
    node.scores <- as.matrix(node.scores)
  }
  if(!is(node.scores, "matrix"))
  {
    if(!is.null(names(node.scores)))
    {
      node.names <- V(network)$name
      node.scores <- node.scores[node.names[which(node.names %in% names(node.scores))]]
    }
    else
    {
      warning("Scores are not named!")
      return(FALSE)
    }
  }
  if(is(node.scores, "matrix"))
  {
    if(!is.null(rownames(node.scores)))
    {
      node.names <- V(network)$name
      node.scores <- node.scores[node.names[which(node.names %in% rownames(node.scores))],]
    }
    else
    {
      warning("Scores are not named!") 
      return(FALSE)
    }
  }
  nodes <- V(network)$name
  # for vector
  if(!is(node.scores, "matrix"))
  {
    table.nodes <- matrix(nrow=length(nodes), ncol=1)
    rownames(table.nodes) <- nodes
    table.nodes[names(node.scores), ] <- node.scores
    table.nodes <- cbind(nodes, table.nodes)
  }
  # for matrix
  if(is(node.scores, "matrix"))
  {
    table.nodes <- matrix(nrow=length(nodes), ncol=dim(node.scores)[2])
    rownames(table.nodes) <- nodes
    table.nodes[rownames(node.scores), ] <- node.scores
    table.nodes <- cbind(nodes, table.nodes)
  }
  c.name=c("#node")
  if(any(is.na(table.nodes)))
  {
    warning("Not all nodes of the network were scored, unscored nodes were set to mean")  
  }
  for(i in 2:dim(table.nodes)[2])
  {
    expect <- sum(as.vector(table.nodes[, i], mode="numeric"), na.rm = TRUE)/length(table.nodes[,i])
    table.nodes[which(is.na(table.nodes[,i])),i] <- expect
    c.name <- c(c.name, paste("score", i-1, sep=""))
  }
  colnames(table.nodes) <- c.name
  write.table(table.nodes, file=paste(file, ".txt", sep=""), append=FALSE, quote=FALSE, sep="\t", dec=".", row.names=FALSE, col.names=TRUE)
  return(TRUE)
}

writeHeinz <- function(network, file, node.scores=0, edge.scores=0, use.node.score=FALSE, use.edge.score=FALSE)
{
  file <- .cleanFile(file)
  if(is(network, "graphNEL"))
  {
    network <- igraph.from.graphNEL(network)
  }
  if(is(node.scores, "data.frame"))
  {
    node.scores <-  as.matrix(node.scores)
  }
  if(is(edge.scores, "data.frame"))
  {
    edge.scores <-  as.matrix(edge.scores)
  }
  # both scores in attributes --> vector
  if(use.node.score == TRUE && use.edge.score==TRUE)
  {
    ret.n <- writeHeinzNodes(network=network, file=paste(file, "_n.txt", sep=""), use.score=use.node.score)
    ret.e <- writeHeinzEdges(network=network, file=paste(file, "_e.txt", sep=""), use.score=use.edge.score)
    return((ret.n&&ret.e))
  }
  # both scores = matrix, check identical dim
  if(is(node.scores, "matrix") && is(edge.scores, "matrix"))
  {
    if(dim(node.scores)[2]<dim(edge.scores)[2])
    {
      warning("Dimension of node.scores smaller than edge.scores, filled with zeros.")
      node.scores <- cbind(node.scores, matrix(0,nrow=dim(node.scores)[1],ncol=dim(edge.scores)[2]-dim(node.scores)[2])) 
    }
    if(dim(node.scores)[2]>dim(edge.scores)[2]) 
    {
      warning("Dimension of edge.scores smaller than node.scores, filled with zeros.")
      edge.scores <- cbind(edge.scores, matrix(0,nrow=dim(edge.scores)[1],ncol=dim(node.scores)[2]-dim(edge.scores)[2])) 
    }
    ret.ne <- writeHeinzNodes(network=network, file=paste(file, "_n.txt", sep=""), node.scores=node.scores)
    ret.n <- writeHeinzEdges(network=network, file=paste(file, "_e.txt", sep=""), edge.scores=edge.scores)
    return((ret.n&&ret.e))
  }
  # one matrix, one vector
  if(is(node.scores, "matrix") || is(edge.scores, "matrix"))
  {
    if(is(node.scores, "matrix"))
    {
      if(use.edge.score)
      {
        edge.scores <- matrix(rep(E(network)$score, dim(node.scores)[2]), nrow=length(E(network)), ncol=dim(node.scores)[2])
      }
      else
      {
        if(length(edge.scores) == 1 && edge.scores == 0)
        {
          edge.scores <- matrix(0, nrow=length(E(network)), ncol=dim(node.scores)[2])
        }
        else
        {
          edge.scores <- matrix(rep(edge.scores, dim(node.scores)[2]), nrow=length(edge.scores), ncol=dim(node.scores)[2])
        }
      }
      ret.n <- writeHeinzNodes(network=network, file=paste(file, "_n.txt", sep=""), node.scores=node.scores, use.score=FALSE)
      ret.e <- writeHeinzEdges(network=network, file=paste(file, "_e.txt", sep=""), edge.scores=edge.scores, use.score=FALSE)
      return(ret.n&&ret.e)
    }
    if(is(edge.scores, "matrix"))
    {
      if(use.node.score)
      {
        node.scores <- matrix(rep(V(network)$score, dim(edge.scores)[2]), nrow=length(V(network)), ncol=dim(edge.scores)[2])
        rownames(node.scores) <- V(network)$name
      }
      else
      {
        if(length(node.scores) == 1 && node.scores == 0)
        {
          node.scores <- matrix(0, nrow=length(V(network)), ncol=dim(edge.scores)[2])
        }
        else
        {
          node.names <- names(node.scores)
          node.scores <- matrix(rep(node.scores, dim(edge.scores)[2]), nrow=length(node.scores), ncol=dim(edge.scores)[2])
          rownames(node.scores) <- node.names
        }
      }
      ret.n <- writeHeinzNodes(network=network, file=paste(file, "_n.txt", sep=""), node.scores=node.scores, use.score=FALSE)
      ret.e <- writeHeinzEdges(network=network, file=paste(file, "_e.txt", sep=""), edge.scores=edge.scores, use.score=FALSE)
      return(ret.n&&ret.e)
    }
  }
  # both scores = vector
  else
  {
    ret.n <- writeHeinzNodes(network=network, file=paste(file, "_n.txt", sep=""), node.scores=node.scores, use.score=use.node.score)
    ret.e <- writeHeinzEdges(network=network, file=paste(file, "_e.txt", sep=""), edge.scores=edge.scores, use.score=use.edge.score) 
    return(ret.n&&ret.e)
  }
  return(FALSE)
}

# run heinz with node and edge files
# run.heinz(heinz.folder="", heinz.e.file, heinz.n.file, N = TRUE, E=FALSE)
# arguments:
#   heinz.folder: folder where Heinz and cplex is stored
#   heinz.e.file: the heinz edge input file, both files are required
#   heinz.n.file: the heinz node input file, both files are required
#   N: boolean, whether to use the node file
#   E: boolean, whether to use the edge file, both files are required
runHeinz <- function (heinz.folder = "", heinz.e.file, heinz.n.file, N = TRUE, 
    E = FALSE, diff = -1, n = 1) 
{
    if (heinz.folder != "") {
        heinz.folder <- file.path(heinz.folder, "", fsep = .Platform$file.sep)
    }
    N = if (N) 
        "True"
    else "False"
    E = if (E) 
        "True"
    else "False"
    command <- paste("cd ", heinz.folder, ";", "heinz.py -e ", heinz.e.file, " -n ", heinz.n.file, " -N ", N, " -E ", E, " -d ", diff, " -s ", n, sep = "")
    system(command)
}
  
#
# *** genereate a dataframe of scores over a given range of fdrs
#
scanFDR <- function(fb, fdr, labels=names(fb$pvalues))
{
  res <- matrix(0, length(fb$pvalues), length(fdr)); 
  colnames(res) <- paste("score_", fdr, sep=""); 
  for(i in 1:length(fdr)) {
    res[,i] <- scoreFunction(fb, fdr[i]);
  }
  if (!is.null(labels)) 
  {
    rownames(res) <- labels 
    return(data.frame(res, stringsAsFactors=F));
  }
  else return(data.frame(res, stringsAsFactors=F));
}


# matrix heinz

# get heinz output, as tree with heinz edges and nodes
# readHeinzTree(heinz.n.output.file, heinz.e.output.file, network)
# arguments:
#   heinz.n.output.file: heinz output file for nodes
#   heinz.e.output.file: heinz output file for edges
#   format: "igraph" or "graphNEL"
# value: igraph object
readHeinzTree <- function(node.file, edge.file, network, format=c("graphNEL", "igraph"))
{
  format <- match.arg(format)
  edges <- as.matrix(read.table(file=edge.file, skip=1, na.strings="n/a"))
  nodes <- as.matrix(read.table(file=node.file, skip=1, na.strings="n/a"))
  modules = list()
  for(i in 2:dim(nodes)[2])
  {
    edges2 <- edges
    nodes2 <- nodes
    if(sum(is.na(as.numeric(edges[,i+1]))) != 0)
    {
      edges2 <- edges[!is.na(as.numeric(edges[,i+1])),]
    }
    if(sum(is.na(as.numeric(nodes[,i]))) != 0)
    {
      nodes2 <- nodes[!is.na(as.numeric(nodes[,i])),]
    }
    if(!is(nodes2, "vector"))
    {
      warning("In at least one result no module found")
      module <- graph.empty(n=0, directed=FALSE)
    }
    else
    {
      if(is(network, "igraph"))
      {
        network <- igraph.to.graphNEL(network)
      }
      # extract complete subgraph
      if(!is(nodes2, "matrix"))
      {
        compl.module <- .subNetwork0(as.character(nodes2[1]), network)
      }
      else
      {
        compl.module <- .subNetwork0(as.character(nodes2[,1]), network)
      }
      rm.edges <- setdiff(edgeNames(compl.module), c(paste(edges2[,1], "~", edges2[,2], sep=""), paste(edges2[,2], "~", edges2[,1], sep="")))
      if(!is.na(rm.edges[1]))
      {
        rm.edges <- matrix(unlist(strsplit(rm.edges, "~")), ncol=length(rm.edges), nrow=2)
        # remove edges that were not found with HEINZ
        module <- (removeEdge(from=rm.edges[1,], to=rm.edges[2,], compl.module))
      }
      else{module <- compl.module}
    }
    if(format == "igraph" && is(module, "graphNEL"))
    {
      module <- igraph.from.graphNEL(module)
    }
    modules[[i-1]] <- module
  }
  if(length(modules)==1)
  {
    return(modules[[1]])
  }
  return(modules)
}

# get heinz output, as graph with heinz nodes and all edges
# readHeinzTree(heinz.n.output.file, heinz.e.output.file, network)
# arguments:
#   heinz.n.output.file: heinz output file for nodes
#   heinz.e.output.file: heinz output fiel for edges
#   network: network that was used for node/edge scores
#   heinz.e.input.file: input edge file that was used for heinz
# value: graph object
readHeinzGraph <- function(node.file, network, format=c("graphNEL", "igraph"))
{
  format <- match.arg(format)
  if(is(network, "graphNEL"))
  {
    network <- igraph.from.graphNEL(network)
  }
  modules <- list()
  # if only node file
  # nodes <- as.matrix(read.table(file=node.file, skip=1, na.strings="n/a"))
  nodesTable <- read.table(file=node.file, skip=1, na.strings="n/a")
  nodesNames <- as.character(nodesTable[,1])
  nodesMatrix <- as.matrix(nodesTable[,seq.int(2, ncol(nodesTable))])
  for(i in seq_len(ncol(nodesMatrix)))
  {
    nodes2 <- nodesMatrix[, i]
    toKeep <- nodesNames[!is.na(as.numeric(nodes2))]
    
    if (length(toKeep) == 0) {
      warning("In at least one result no module found")
      module <- graph.empty(n=0, directed=FALSE)
    }
    else
    {
      module <- .subNetwork0(toKeep, network)
    }
    if(format == "graphNEL")
    {
      module <- igraph.to.graphNEL(module)
    } 
    modules[[i]] <- module
  }
  if(length(modules)==1)
  {
    return(modules[[1]])
  }
  return(modules)      
}



# heuristic to get max. scoring subnetwork
runFastHeinz <- function(network, scores) 
{
    net.flag <- FALSE
    if(is.null(names(scores))) {
        warning("Unnamed scores")
    }
    if(any(is.na(names(scores)))){
        warning("NA in names of scores, score will be removed")
        scores <- scores[!is.na(names(scores))]
    }
    if (is(network, "graphNEL")) {
        network <- igraph.from.graphNEL(network)
        net.flag <- TRUE
    }
    if(is.null(V(network)$name))
    {
      V(network)$name <- as.character(V(network))
    }
    V(network)$score <- scores[V(network)$name]
    pos.nodes <- names(scores[which(scores > 0)])
    if (length(pos.nodes) == 0) {
        warning("No positive nodes")
        module <- graph.empty(n = 0, directed = FALSE)
        if (net.flag) {
			nE <- ecount(module)
			module <- simplify(module, remove.multiple = TRUE)
			if(nE!=ecount(module))
			{
				warning("Multiple edges between two nodes had to be removed for calculation")
			}
            module <- igraph.to.graphNEL(module)
        }
        return(module)
    }
    if (length(pos.nodes) == 1) {
        module <- .subNetwork0(pos.nodes, network)
        if (net.flag) {
			nE <- ecount(module)
			module <- simplify(module, remove.multiple = TRUE)
			if(nE!=ecount(module))
			{
				warning("Multiple edges between two nodes had to be removed for calculation")
			}
            module <- igraph.to.graphNEL(module)
        }
        return(module)
    }
    pos.subgraph <- .subNetwork0(pos.nodes, network)
    conn.comp.graph <- decompose.graph(pos.subgraph)
    score.comp <- unlist(lapply(lapply(conn.comp.graph, get.vertex.attribute, "score"), sum))
    conn.comp.graph <- conn.comp.graph[order(score.comp, decreasing = TRUE)]
    score.comp <- sort(score.comp, TRUE)
    for (i in 1:length(conn.comp.graph)) {
        conn.comp.graph[[i]]$score <- score.comp[i]
    }
    v.id <- seq(1, vcount(network))
    names(v.id) <- V(network)$name
    edgelist <- get.edgelist(network, FALSE)
    edgelist1 <- edgelist[, 1]
    edgelist2 <- edgelist[, 2]
    for (i in 1:length(conn.comp.graph)) {
        new.id <- length(V(network)) + i
        for (j in as.character(v.id[V(conn.comp.graph[[i]])$name])) {
            edgelist1[which(edgelist1 == j)] <- new.id
            edgelist2[which(edgelist2 == j)] <- new.id
        }
    }
    new.ids <- seq(length(V(network)) + 1, length(V(network)) + length(conn.comp.graph))
    new.names <- paste("cluster", seq(1:length(conn.comp.graph)), sep = "")
    names(new.ids) <- new.names
    v.id <- c(v.id, new.ids)
    v.name <- names(v.id)
    names(v.name) <- v.id
    new.edgelist <- cbind(v.name[as.character(edgelist1)], v.name[as.character(edgelist2)])
    interactome2 <- graph.edgelist(new.edgelist, FALSE)
    E(interactome2)$weight <- rep(0, length(E(interactome2)))
    interactome2 <- simplify(interactome2, remove.loops = TRUE, remove.multiple = TRUE)
    score1 <- scores[V(interactome2)$name]
    names(score1) <- V(interactome2)$name
    score1[which(is.na(score1))] <- 0
    score.degree <- score1/(degree(interactome2) + 1)
    V(interactome2)$score.degree <- score.degree
    E(interactome2)$weight <- -(V(interactome2)[get.edgelist(interactome2, FALSE)[, 1]]$score.degree + V(interactome2)[get.edgelist(interactome2, FALSE)[, 2]]$score.degree)
    node.score <- scores[V(interactome2)$name]
    names(node.score) <- V(interactome2)$name
    node.score.cluster <- sapply(conn.comp.graph, get.graph.attribute, "score")
    names(node.score.cluster) <- new.names
    node.score[grep("cluster", names(node.score))] <- node.score.cluster[names(node.score[grep("cluster", names(node.score))])]
    decomp.graphs <- decompose.graph(interactome2)
    sum.pos <- lapply(decomp.graphs, function(x) {
        sum(node.score[names(which(node.score[V(x)$name] > 0))])
    })
    interactome2 <- decomp.graphs[[which.max(sum.pos)]]
    rm(decomp.graphs)
    mst <- minimum.spanning.tree(interactome2, weights = E(interactome2)$weight)
    mst.cluster.id <- grep("cluster", V(mst)$name)
    names(mst.cluster.id) <- V(mst)[mst.cluster.id]$name
    mst.cluster.id <- mst.cluster.id[order(as.numeric(matrix(unlist(strsplit(names(mst.cluster.id), "cluster")), nrow = 2)[2, ]))]
    all.ids <- c()
		# check if multiple clusters exist, 
		# if they do not exist, assign neg nodes to empty vector that is used later.
		if(length(mst.cluster.id)==1){
			neg.node.ids.2 = c()
		}else
		{
			for (j in 1:(length(mst.cluster.id)-1)) {
					path <- unlist(get.all.shortest.paths(mst, from = mst.cluster.id[j], to = mst.cluster.id[(j + 1):length(mst.cluster.id)]))
					all.ids <- c(all.ids, path)
			}
			all.ids <- unique(all.ids)
			sub.interactome2 <- .subNetwork0(V(mst)[all.ids]$name, interactome2)
			neg.node.ids <- which(node.score[V(sub.interactome2)$name] < 0)
			for (i in neg.node.ids) {
					V(sub.interactome2)[i]$clusters <- list(neighbors(sub.interactome2, v = i)[grep("cluster", V(sub.interactome2)[neighbors(sub.interactome2, v = i)]$name)])
			}
			score.neg.nodes <- c()
			for (i in neg.node.ids) {
			    # >>> FIX !!! Return values were changed from NA to empty vertex set
			    #if(!is.na(V(sub.interactome2)[i]$clusters[[1]][1]))
			    # >>> FIX !!!
			    if (length(V(sub.interactome2)[i]$clusters[[1]])>0) {
						score.neg.nodes <- c(score.neg.nodes, sum(node.score[V(sub.interactome2)[c(i, V(sub.interactome2)[i]$clusters[[1]])]$name]))
					}
					else {
						score.neg.nodes <- c(score.neg.nodes, node.score[V(sub.interactome2)[i]$name])
					}
			}
			neg.node.ids.2 <- neg.node.ids[score.neg.nodes > 0]
		}
    if (length(neg.node.ids.2) == 0) {
        module <- unlist(lapply(conn.comp.graph, get.vertex.attribute, "name")[as.numeric(matrix(unlist(strsplit(names(node.score.cluster)[which.max(node.score.cluster)], "cluster")), nrow = 2)[2,])])
        module <- .subNetwork0(module, network)
        if (net.flag) {
			nE <- ecount(module)
			module <- simplify(module, remove.multiple = TRUE)
			if(nE!=ecount(module))
			{
				warning("Multiple edges between two nodes had to be removed for calculation")
			}
            module <- igraph.to.graphNEL(module)
        }
        return(module)
    }
    subg <- largestComp(induced.subgraph(sub.interactome2, neg.node.ids.2))
    mst.subg <- minimum.spanning.tree(subg, E(subg)$weight)
    max.score <- 0
    best.path <- c()
    for (i in 1:(length(V(mst.subg)))) {
        path <- get.all.shortest.paths(mst.subg, from = V(mst.subg)[i])
        path.score <- unlist(lapply(path$res, .getPathScore, graph1 = mst.subg, graph2 = sub.interactome2, node.score = node.score))
        best.pos <- which.max(path.score)
        if (path.score[[best.pos]] > max.score) {
            best.path <- path$res[[best.pos]]
            max.score <- path.score[[best.pos]]
        }
    }
    if(length(best.path)!=1){
      cluster.list <- V(mst.subg)[best.path]$clusters
  	  names.list <- as.character(1:length(cluster.list))
  	  names(cluster.list) <- names.list
  	  names(best.path) <- names.list

      for (i in names.list) {
          res <- lapply(cluster.list, intersect, cluster.list[[i]])
  		  if(length(intersect(unlist(cluster.list[as.character(which(as.numeric(names.list)<as.numeric(i)))]), unlist(cluster.list[as.character(which(as.numeric(names.list)>as.numeric(i)))])))>0){
  			  if (length(setdiff(res[[i]], unique(unlist(res[names(res)!=i]))))==0){
  				  cluster.list <- cluster.list[names(cluster.list)!=i]
  				  names.list <- names.list[names.list!=i]
  			  }
        }
      }
      best.path <- best.path[names.list]
    }
    module <- V(mst.subg)[best.path]$name
    pos.cluster <- V(sub.interactome2)[unique(unlist(V(mst.subg)[best.path]$clusters))]$name
    module <- c(module, unlist(lapply(conn.comp.graph, get.vertex.attribute, "name")[as.numeric(matrix(unlist(strsplit(pos.cluster, "cluster")), nrow = 2)[2, ])]))
    module <- .subNetwork0(module, network)
    if (net.flag) {
		nE <- ecount(module)
		module <- simplify(module, remove.multiple = TRUE)
		if(nE!=ecount(module))
		{
			warning("Multiple edges between two nodes had to be removed for calculation")
		}
        module <- igraph.to.graphNEL(module)
    }
    return(module)
}


# function for lapply
.getPathScore <- function(path, graph1, graph2, node.score)
{
  sum(c(node.score[V(graph1)[path]$name], node.score[V(graph2)[unique(unlist(V(graph1)[path]$clusters))]$name]))
}

