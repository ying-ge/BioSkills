
# *********************************************************
# *
# * Statistics: Order statistics and BUM Model
# *
# *
#
#
#
#
#
# *********************************************************

# aggr.pvals(pval.array, order, plot=TRUE)
# arguments: 
#   pval.array: array of p-values for several experiments/covariates
#   order: the order statistic which is used
#   plot: boolean, whether to plot the aggregated p-value distributions
# value: vector of aggregated p-values of the given order
aggrPvals <- function(pval.matrix, order=ncol(pval.matrix), plot=TRUE)
{
  old.pvals <- pval.matrix
  if(!is.matrix(pval.matrix))
  {
    return("Input is not matrix")
  }
  if(order>dim(pval.matrix)[2])
  {
    return("order is larger than array dimensions")
  }
  for(j in 1:dim(pval.matrix)[1])
  {
    pval.matrix[j,] <- sort(as.numeric(pval.matrix[j,]), na.last=TRUE)
  }
  x.vec <- as.numeric(pval.matrix[,order])
  n <- dim(pval.matrix)[2]
  concat.pvals <- pbeta(x.vec, order, n-order+1)
  names(concat.pvals)<- rownames(pval.matrix)
  if(plot)
  {
    for(j in 1:dim(pval.matrix)[2])
    {
      hist(as.numeric(old.pvals[,j]), n=50, main=paste("Histogram of ", j, ". p-values", sep=""))
      par(ask = TRUE)
    }
    hist(as.numeric(concat.pvals), n=50, main="Histogram of aggregated p-values")
  }
  return(concat.pvals)
}

# 
# *** beta uniform mixture model with shape2 fixed @1
#
fbum <- function(x, lambda, a){ lambda+(1-lambda)*a*x^(a-1) };

#
#     Log likelihood of BUM model
# ... version for optim
fbumLL  <- function(parms, x){sum(log(fbum(x, parms[1], parms[2])))};

# ... standard optim
bumOptim <- function(x, starts=1, labels=NULL)
{
  if(is.null(names(x)) && is.null(labels))
  {
    warning("Please name the p-values with the gene names or give labels!")
    names(x) <- as.character(1:length(x))
  }
  if(!is.null(labels))
  {
	names(x) <- labels
  }
  a <- runif(starts, 0.3, 0.7)
  lambda <- runif(starts, 0.3, 0.7)
  value <- Inf
  best <- list()
  for(i in 1:starts)
  {
    test.optim <- try(opt <- optim(c(lambda[i], a[i]), fn=.fbumnLL, gr=.fpLL, x=x, lower=rep(1e-5,3), method="L-BFGS-B", upper=rep(1-1e-5,3)))
    if ((!class(test.optim)=="try-error") && all(opt$par >= 1e-5) && all(opt$par <= 1-1e-5))
    {
      value <- opt$value
      best <- opt
    }
  }
  if(length(best)==0)
  {
    return(warning("BUM model could not be fitted to data"))
  }
  else
  {
    if (any(opt$par == 1e-5) || any(opt$par == 1-1e-5))
    {
      warning("One or both parameters are on the limit of the defined parameter space")
    }
    ret <- list(lambda=best$par[1], a=best$par[2], negLL=best$value, pvalues=x)
    class(ret) <- "bum"
    return(ret)
  }
}

# print function
print.bum <- function(x, ...)
{
  cat("Beta-Uniform-Mixture (BUM) model\n\n");
  cat(paste(length(x$pvalues), "pvalues fitted\n\n"));
  cat(sprintf("Mixture parameter (lambda):\t%1.3f\n", x$lambda));
  cat(sprintf("shape parameter (a): \t\t%1.3f\n", x$a));
  cat(sprintf("log-likelihood:\t\t\t%.1f\n", -x$negLL));
}

# summary function
summary.bum <- function(object, ...)
{
  cat("Beta-Uniform-Mixture (BUM) model\n\n");
  cat(paste(length(object$pvalues), "pvalues fitted\n\n"));
  cat(sprintf("Mixture parameter (lambda):\t%1.3f\n", object$lambda));
  cat(sprintf("shape parameter (a): \t\t%1.3f\n", object$a));
  cat(sprintf("log-likelihood:\t\t\t%.1f\n", -object$negLL));
  cat(sprintf("pi-upper bound:\t\t\t%.1f\n", piUpper(object)));
}

# fit bum model to p-values
# fit.fb(p.values, plot=TRUE)
# arguments:
#   p.values: p-values
#   plot: whether to plot a qqplot and a histogram of the fitted values
# values: fitted model
fitBumModel <- function (x, plot = TRUE, starts=10) 
{
    if (is.null(names(x)))
    {
        warning("Please name the p-values with the gene names!")
        names(x) = as.character(1:length(x))
    }
    fit <- bumOptim(x = x, starts)
    if (plot) 
    {
        par(mfrow = c(1, 2))
        hist(x = fit)
        plot(fit)
    }
    return(fit)
}

#
# calculate tau from given FDR
#
fdrThreshold <- function(fdr, fb)
{
  pihat <- fb$lambda+(1-fb$lambda)*fb$a;
  return(((pihat - (fdr*fb$lambda)) / (fdr*(1-fb$lambda)))^(1/(fb$a-1)))
}

#
# *** upper bound pi for the fraction of noise
#     (see Pounds an Morris)
#
piUpper <- function(fb) { return(fb$lambda + (1 - fb$lambda)*fb$a) }


# ... offset constant for score
#     changing from FDR1 to FDR2
scoreOffset <- function(fb, fdr1, fdr2)
{
  return((fb$a-1)* (log(fdrThreshold(fdr1, fb))- log(fdrThreshold(fdr2, fb))))
}

#
#
#
scoreFunction <- function(fb, fdr=0.01)
{
  return((fb$a - 1)*(log(fb$pvalues)-log(fdrThreshold(fdr, fb))))
}

#
# *** Score the nodes of a network 
#     
#     TODO: fill score into score attributes 
#
scoreNodes <- function(network, fb, fdr=0.05)
{
  score <- scoreFunction(fb, fdr);
  if(is(network, "graphNEL"))
  {
    score <- score[match(nodes(network),names(fb$pvalues))];
    names(score) <- nodes(network);
    return(score);
  }
  else if(is(network, "igraph"))
  {
	if(is.null(V(network)$name))
    {
      V(network)$name <- as.character(V(network))
    }
    score <- score[V(network)$name]
    return(score)
  }
}


#
# *** returns a dataframe with scores and components
#
getCompScores <- function(network, score)
{
  if(is(network, "igraph"))
  {
    cl.network <- clusters(network) 
    cc <- order(cl.network$csize, decreasing=TRUE)
    cc <- lapply(cc, memb= cl.network$membership, function(x, memb=memb) V(network)[which(memb==x)]$name)
  }
  else
  {
	 cc <- connComp(network)   
	 cc <- cc[order(listLen(cc), decreasing=TRUE)];  
  }
  len <- listLen(cc)
  label <- unlist(cc);
  comp <- unlist(mapply(rep,1:length(len), len ));
  score <- score[match(label, names(score))];
  return(data.frame(comp=comp,label=label, score=score));
}


# *********************************************************
# *
# * internal  
# *
# *********************************************************


#
#
#
# ... gradient of .fLL
.fpLL  <-function(parms, x)
{
 l <- parms[1]; a <- parms[2];
 
 dl <- -sum((1-a*x^(a-1))/(a*(1-l)*x^(a-1)+l));
 
 da <- -sum((a*(1-l)*x^(a-1)*log(x)+(1-l)*x^(a-1))/(a*(1-l)*x^(a-1)+l));

 return(c(dl,da));
}

# negative log likelihood
.fbumnLL <- function(parms, x){-fbumLL(parms, x)}



