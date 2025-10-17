
# *********************************************************
# *
# * Plots
# *
# *********************************************************

# qqplot.bum(pvalues, fb)
# arguments:
#   pvalues: vector of p-values
#   fb: fitted bum model to the p-value distribution
# values: plot -> quantiles of the bum distribution and the observed p-values
plot.bum <- function(x, main="QQ-Plot", xlab="Estimated p-value", ylab="Observed p-value", ...)
{
	n <- length(x$pvalues)
	probs <- (rank(sort(x$pvalues))-.5)/n
	# get quantiles of the bum distribution
	quantiles <- unlist(sapply(probs, uniroot, f=.pbum.solve, interval=c(0,1), lambda=x$lambda, a=x$a)[1,])
	plot(c(0,1),c(0,1), main=main, xlab=xlab, ylab=ylab, type="n", ...)
	lines(quantiles, sort(x$pvalues), lty=2) 
	lines(c(0,1),c(0,1), col="grey")		
}

# hist.bum(pvalues, fb)
# arguments:
#   pvalues: vector of p-values
#   fb: fitted bum model to the p-value distribution
# values: plot -> histogram of p-values with fitted bum distribution
hist.bum <- function(x, breaks=50, main="Histogram of p-values", xlab="P-values", ylab="Density", ...)
{
	hist(x$pvalues, breaks=breaks,  probability=TRUE, main=main, xlab=xlab, ylab=ylab, ...)
	bum.data <- seq(from=0, to=1, 1/100)
	lines(bum.data, x$lambda+(1-x$lambda)*x$a*bum.data^(x$a-1), lwd=3, col="red3");
	abline(h=piUpper(x), col="blue3", lwd=2);
	axis(side=2, labels=expression(pi), at=piUpper(x));
}

# internal function for root detection
.pbum.solve <- function(x, lambda ,a, proba)
{
  return((lambda*x+(1-lambda)*x^a)-proba)
}


# *** nice graphics...
plotLLSurface <- function(x, opt=NULL, main="Log-Likelihood Surface", color.palette=heat.colors, nlevels=32)
{
  if(is.null(opt)) opt <- list(a=-1.0, l=-1.0); 
  f   <- function(l, a) {fbumLL(c(l, a), x)};
  
  v <- seq(0.1, 0.9, 0.05);
  z <- outer(v, v, Vectorize(f));
  
  Lines <- list(bquote(lambda == .(round(opt$l, 4))), bquote(a == .(round(opt$a, 4))))
  
  filled.contour(v, v, z, nlevels=nlevels, color.palette=color.palette,
  main=main, xlab=expression(lambda), ylab="a",
  plot.axes={axis(1, seq(0,1,0.1)); axis(2, seq(0,1,0.1));        
    abline(v=opt$l, lty=2, col="darkgray");
    abline(h=opt$a, lty=2, col="darkgray");
    hght <- strheight("Here")
    points(opt$l, opt$a, cex=1.5);
    text(opt$l, opt$a - (hght * 1.5*seq(length(Lines))), do.call(expression, Lines), adj=c(-0.2, -2.3), cex=c(0.8, 0.8));
    }); 
} 


# graph plot in igraph
plotModule <- function (network, layout = layout.fruchterman.reingold, labels = NULL, diff.expr = NULL, scores = NULL, main = NULL, vertex.size=NULL, ...) 
{    
    if (is(network, "graphNEL"))
    {
        network <- igraph.from.graphNEL(network)
    }
    if (is.null(V(network)$name))
    {
        V(network)$name <- as.character(V(network))
    }
    if (is.null(labels))
    {
        if ("geneSymbol" %in% list.vertex.attributes(network))
        {
            labels <- V(network)$geneSymbol
        }
        else
        {
            labels <- V(network)$name
        }

    }
    shapes <- rep("circle", length(V(network)))
    names(shapes) <- V(network)$name
    if (!is.null(scores) && !is.null(names(scores)))
    {
        shapes[intersect(names(which(scores < 0)), V(network)$name)] <- "csquare"
    }
    if (is.null(scores) && "score" %in% list.vertex.attributes(network))
    {
        scores <- V(network)$score
        names(scores) <- V(network)$name
        shapes[names(which(scores < 0))] <- "csquare"
    }
    if (!is.null(diff.expr) && !is.null(names(diff.expr)))
    {
        coloring <- .node.color(network, diff.expr)
    }
    else
    {
        coloring <- "SkyBlue2"
    }
    if (is.null(diff.expr) && "diff.expr" %in% list.vertex.attributes(network))
    {
        diff.exprs = V(network)$diff.expr
        names(diff.exprs) <- V(network)$name
        coloring <- .node.color(network, diff.exprs)
    }
    max.labels <- max(nchar(labels))
    network.size = length(V(network))
    vertex.size2 <- 8
    cex = 0.6
    if (network.size < 50)
    {
        if (max.labels > 2)
        {

            labels.dist <- 0.5
        }
        else
        {
            vertex.size2 <- 15
            labels.dist <- 0
        }
    }
    if (network.size < 100 && network.size >= 50)
    {

        if (max.labels > 2)
        {
            labels.dist <- 0.5
        }
        else
        {
            labels.dist <- 0
        }
    }
    if (network.size >= 100)
    {

        if (max.labels > 3)
        {
            labels.dist <- 0.5
        }
        else
        {
            labels.dist <- 0
        }
    }
    if(!is.null(vertex.size))
    {
      vertex.size2 <- vertex.size
      labels.dist <- vertex.size/15
    }
    plot(network, layout = layout, vertex.size = vertex.size2, 
        vertex.label = labels, vertex.label.cex = cex, vertex.label.dist = labels.dist, 
        vertex.color = coloring, vertex.label.family = "sans", 
        vertex.shape = shapes, main = main, ...)
}


# Plot the network in 3D
plot3dModule <- function(network, labels=NULL, windowSize = c(100,100,1500,1000), diff.or.scores=NULL, red=c("negative", "positive"), ...)
{
  rgl.loaded <- requireNamespace("rgl")
  if(!rgl.loaded){warnings("Please install rgl package for this plot")}
  else
  {
    red <- match.arg(red)
    if(is(network, "graphNEL"))
    {
      network <- igraph.from.graphNEL(network)
    }
    if(is.null(V(network)$name))
    {
      V(network)$name <- seq(from=1, to=length(V(network)))
    }
    if(is.null(labels))
    {
      if("geneSymbol" %in% list.vertex.attributes(network))
      {
        labels <- V(network)$geneSymbol
      }
      else{ labels <- V(network)$name}
    }
    if(!is.null(diff.or.scores) && !is.null(names(diff.or.scores)))
    {
      if(red == "negative")
      {
        diff.or.scores <- -diff.or.scores
      }
      coloring <- .node.color(network, diff.or.scores)
    }  
    else
    { 
      coloring <- "SkyBlue2"
    }
    if(is.null(diff.or.scores) && "score" %in% list.vertex.attributes(network))
    {
      scores = V(network)$score
      names(scores) <- V(network)$name
      if(red == "negative")
      {
        scores <- -scores
      }
      coloring <- .node.color(network, scores)
    } 
    max.labels <-  max(nchar(labels))
    network.size = length(V(network))
    coords <- layout.fruchterman.reingold(network, dim=3) 
    if(network.size<50)
    {
      vertex.size <- 10
      if(max.labels > 3)
      {
        labels.dist <- 0.5
      } 
      else{labels.dist <- 1}
    }
    if(network.size<100 && network.size>=50)
    {
      vertex.size <- 8
      if(max.labels > 3)
      {
        labels.dist <- 0.3
      } 
      else{labels.dist <- 0.5}
    }
    if(network.size>=100)
    {
      vertex.size <- 6
      if(max.labels > 3)
      {
        labels.dist <- 0.3
      } 
      else{labels.dist <- 0.5}
    }  
    rgl.open()
    par3d(windowRect=windowSize, family="sans", zoom=0.7)
    rgl.texts(x=c(0,0,0), text="")
    rglplot(network, layout=coords, vertex.size=vertex.size, vertex.color=coloring, vertex.label=labels, vertex.label.dist=labels.dist, vertex.label.family="sans", ...)
    rgl.bg(sphere=TRUE, color="lightskyblue", back="filled")
  }
}

save3dModule <- function(file)
{
  requireNamespace("rgl")
  file <- .cleanFile(file)
  rgl.bg(color="white")
  rgl.postscript(filename=paste(file, ".pdf", sep=""), fmt="pdf")
}

# internal methods #############################################################

.node.color <- function(network, colors)
{  
  colors <- colors[V(network)$name]
  colors2 <- colors
  # set red colors
  if(max(abs(colors))<5)
  {
    colors <- colors*5
  }
  if(any(colors>0))
  {
	  max.red <- max(ceiling(abs(colors[which(colors>0)])))
	  reds <- colorRampPalette(colors=c("white", "red"))
	  red.vec <- reds(max.red)
	  colors2[which(colors>0)] <- red.vec[ceiling(abs(colors[which(colors>0)]))]
	  }
  # set green colors
  if(any(colors<0))
  {
	  max.green <- max(ceiling(abs(colors[which(colors<0)])))
	  greens <- colorRampPalette(colors=c("white", "green"))
	  green.vec <- greens(max.green)
	  colors2[which(colors<0)] <- green.vec[ceiling(abs(colors[which(colors<0)]))]
  }
  return(colors2)
}

