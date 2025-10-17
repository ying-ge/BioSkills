
# *********************************************************
# *
# * MAUtils: Utilities for microarray analysis
# *
# *********************************************************

# .affyID2ppiID(affyID)
# arguments:
#   affy ID
# value: human PPI ID
.affyID2ppiID <- function(exprSet, network=NULL, attr="geneID")
{  
  package <- paste(annotation(exprSet), ".db", sep="")
  # check if package is allready loaded
  package.loaded <- paste("package:", package, sep="") %in% search()
  loaded.here <- FALSE
  if(!package.loaded)
  {
    loaded.here <- require(package=package, character.only=TRUE) 
    # check if package was loaded, else return FALSE and warning
    if(!loaded.here)
    {
      warning(paste("Please install ", annotation(exprSet), ".db first", sep=""))
      return(NULL)
    }
  }
  # get annotations for chip
  id <- AnnotationDbi::mget(featureNames(exprSet), get(paste(annotation(exprSet), "ENTREZID", sep="")))
  if(is.null(network))
  {
    symbol <- AnnotationDbi::mget(featureNames(exprSet), get(paste(annotation(exprSet), "SYMBOL", sep="")))
    ppiID <- paste(symbol, "(", id, ")", sep="")
    names(ppiID) <- featureNames(exprSet)
  }
  else
  {
    if(is(network, "igraph"))
    {
	  if(is.null(V(network)$name))
	  {
        V(network)$name <- as.character(V(network))
	  }
      nodes <- V(network)$name
      names(nodes) <- get.vertex.attribute(network, name=attr, index=V(network))
    }
    if(is(network, "graphNEL"))
    {
      nodes <- nodes(network)
      names(nodes) <- unlist(nodeData(network, n=nodes(network), attr=attr))
    }
    # map over entrez IDs, because symbols are not always the same
    ppiID <- nodes[unlist(id)]   
    names(ppiID) <- featureNames(exprSet) 
  }
  # detach package if it was not already loaded
  if(!package.loaded){detach(pos=which(search() == paste("package:", package, sep="")))}
  return(ppiID) 
}

# mapByVar(exprSet)
# arguments:
#   exprSet: expressionSet
#   ignoreAFFX: TRUE / FALSE
# value: returns matrix with one gene (ppiID) per probeset, with the highest variance in the intensity
mapByVar <- function(exprSet, network=NULL, attr="geneID", ignoreAFFX=TRUE)
{
  intens.matrix <- exprs(exprSet)
  variance <- apply(intens.matrix, 1, var)
  var <- sort(variance, TRUE)
  sorted.data <- intens.matrix[names(var),]

  # get ppiIDs for affyIDs
  ppi.names <- .affyID2ppiID(exprSet, network, attr)
  if(is.null(ppi.names))
  {
    return(NULL)
  }
  # remove probes, that could not be mapped (with network)
  pos <- which(is.na(ppi.names))
  if(length(pos)>0)
  {
    ppi.names <- ppi.names[-pos]
  }
  # remove probes, that could not be mapped (without network)
  pos <- which(ppi.names == "NA(NA)")
  if(length(pos)>0)
  {
    ppi.names <- ppi.names[-pos]
  }
  # remove AFFX genes
  if(ignoreAFFX)
  {
    pos <- grep("AFFX", names(ppi.names))
    if(length(pos)>0)
    {
      ppi.names <- ppi.names[-pos]
    }
  }
  sorted.data2 <- sorted.data[names(ppi.names),]
  rownames(sorted.data2) <- ppi.names

  # get first occurence of ppiID
  sorted.data2 <- sorted.data2[match(unique(ppi.names), rownames(sorted.data2)),]
  return(sorted.data2)
}