
# *********************************************************
# *
# * Resampling functions
# *
# *********************************************************

# resamplingPvalues(exprMat, groups, alternative = c("two.sided", "less", "greater"), resampleMat=FALSE)
# arguments:
#   exprMat: expression matrix
#   groups: groups for t.test
#   alternative: alternative for t.test
#   resampeMat: Boolean, return expression matrix with resampled values
# values: list of p-values and resampleMat
resamplingPvalues <- function(exprMat, groups, alternative = c("two.sided", "less", "greater"), resampleMat=FALSE)
{
  alternative <- match.arg(alternative)
  if(length(levels(groups))==2)
  {
    group1 <- which((groups)==levels(groups)[1])
    group2 <- which((groups)==levels(groups)[2])
    resMat <- c()
    if(length(group1)>4 && length(group2)>4)
    {
      jack.p.values <- c()
      size1 <- ceiling(length(group1)*0.5)
      size2 <- ceiling(length(group2)*0.5)
      jack.mat1 <- t(apply(exprMat, 1, function(x) x[sample(group1, size=size1, replace=FALSE)]))
      jack.mat2 <- t(apply(exprMat, 1, function(x) x[sample(group2, size=size2, replace=FALSE)]))
      if(alternative=="two.sided")
      {
        requireNamespace("genefilter")
        jack.p.values <- rowttests(cbind(jack.mat1,jack.mat2), fac=factor(c(rep("A", size1), rep("B", size2))))$p.value
	names(jack.p.values) <- rownames(exprMat)
      }
      else
      {
	      n <- size1 + size2
	      mean1 <- rowMeans(jack.mat1, na.rm=TRUE)
	      mean2 <- rowMeans(jack.mat2, na.rm=TRUE)
	      var1 <- apply(jack.mat1, 1, var, na.rm=TRUE)
	      var2 <- apply(jack.mat2, 1, var, na.rm=TRUE)
	      vars <- (size1-1)/(n-2)*var1 + (size2-1)/(n-2)*var2
	      tstat <- sqrt(size1*size2/n)*(mean1-mean2)/sqrt(vars)
	      jack.p.values <- switch(alternative,
			#two.sided = pt(abs(tstat), df=(n-1), lower.tail = FALSE)*2,
			less = pt(tstat, df=(n-1), lower.tail = TRUE),
			greater = pt(tstat, df=(n-1), lower.tail = FALSE))
      }
      if(resampleMat)
      {
        resMat <- cbind(jack.mat1, jack.mat2)
      }
      return(list(p.values = jack.p.values, resampleMat=resMat))
    }
    else stop("Too few sampled per group, at least 4 samples are needed for each group")
  }
  else stop("Only for two-group camparisons")
}

# consensusScores(modules, network, ro=length(modules)/2)
# arguments:
#   modules: modules from pseudo-replicates of expression values
#   network: interaction network, which should be scores
#   ro: threshold which is subtracted from the scores
# values: list of node scores, edge scores, node frequencies and edge frequencies
consensusScores <- function(modules, network, ro=length(modules)/2)
{
  if(is(modules[[1]], "graphNEL") && is(network, "graphNEL"))
  {
    node.counts <- table(unlist(lapply(modules, function(x) nodes(x))))
    node.scores <- rep(0, numNodes(network))
    names(node.scores) <- nodes(network)
    node.scores[names(node.counts)] <- node.counts
    node.scores2 <- node.scores-ro
    edge.counts <- table(unlist(lapply(lapply(modules, function(x) x), function(x) apply(getEdgeList(x)[, c("from", "to")], 1,  function(y) paste(sort(y), collapse="|")))))
    edge.scores <- rep(0, numEdges(network))
    names(edge.scores) <- apply(getEdgeList(network)[, c("from", "to")], 1,  function(y) paste(sort(y), collapse="|"))
    edge.scores[names(edge.counts)] <- edge.counts
    edge.scores2 <- edge.scores-ro
    return(list(N.scores=node.scores2, E.scores=edge.scores2, N.freq=node.scores/length(modules), E.freq=edge.scores/length(modules)))
  }
  else if(is(modules[[1]], "igraph") && is(network, "igraph"))
  {
    node.counts <- table(unlist(lapply(modules, function(x) V(x)$name)))
    node.scores <- rep(0, vcount(network))
    names(node.scores) <- V(network)$name
    node.scores[names(node.counts)] <- node.counts
    node.scores2 <- node.scores - ro
    edge.counts <- table(unlist(lapply(lapply(modules, function(x) x), function(x) apply(get.edgelist(x), 1,  function(y) paste(sort(y), collapse="|")))))
    edge.scores <- rep(0, ecount(network))
    names(edge.scores) <- apply(get.edgelist(network), 1,  function(y) paste(sort(y), collapse="|"))
    edge.scores[names(edge.counts)] <- edge.counts
    edge.scores2 <- edge.scores - ro
    return(list(N.scores=node.scores2, E.scores=edge.scores2, N.freq=node.scores/length(modules), E.freq=edge.scores/length(modules)))
  }
  else stop("Please give modules and network in same graph format")
}

# sortedEdgeList(network)
# arguments:
#   network: an interaction network
# values: returns a sorted edgelist where the source protein is alphabetically smaller than the target protein
sortedEdgeList <- function(network)
{
  if(is(network, "graphNEL"))
  {
    if(isDirected(network))
    {
      stop("Method only available for undirected networks")
    }
    return(apply(getEdgeList(network)[, c("from", "to")], 1,  function(y) paste(sort(y), collapse="|")))
  }
  if(is(network, "igraph"))
  {
    if(is.directed(network))
    {
      stop("Method only available for undirected networks")
    }
    return(apply(get.edgelist(network), 1,  function(y) paste(sort(y), collapse="|")))
  }
}



