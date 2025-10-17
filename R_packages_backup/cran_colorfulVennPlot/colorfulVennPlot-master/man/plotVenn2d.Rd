\name{plotVenn2d}
\alias{plotVenn2d}
\title{Plot Venn diagram for 2-dimensional data.}
\description{Given a vector of 3 values, which describe 2-dimensional data, it plots a Venn diagram, i.e. 'crossing circles'. The user can specify values, labels for each circle-group and colors.
}
\usage{
plotVenn2d(x, labels = c('A', 'B'),
  Colors = c("red", "yellow", "green"),
  Title = NULL, shrink = 1, rot=0, radius= c(1,1), resizePlot = 1, 
  reverseLabelOrdering=TRUE)
}
\arguments{
  \item{x}{ a numeric vector of length 3, where first value corresponds to only group 1, second is only group 2, and third is a crossover of the two groups. Can also be a character vector of length 3.}
  \item{labels}{ a character vector of length 2, providing names for the 2 groups/dimensions. }
  \item{shrink}{ a numeric value, specifying zooming effect of the plot, defaults to 1. }
  \item{Colors}{ a vector of color names for the backgrounds of each part of the diagram. }
  \item{Title}{ optional: a character vector of length 1, specifying title for the whole plot. }
  \item{rot}{ a numeric value, specifying the number of degrees to rotate the graph. }
  \item{radius}{ a 2-item numeric vector containing the relative sizes of the two circles }
  \item{resizePlot}{ a numeric value indicating the amount to increase or decrease the size of the plot }
  \item{reverseLabelOrdering}{ boolean value indicating whether labels are in the normal or reverse ordering. default=TRUE for backward compatability with the original version. }
}
\details{
NOTICE: This only works for 2-dimensional data.

By definition, value '00' lies outside the plotted diagram. Therefore, it will not be plotted. Because of that, if plotting 'percent' data, all values will not add up to 100 percent, missing 00's allocation.
}
\author{Elliot Noma, Aliona Manvae
}
\note{ It requires package grid to be installed, and can be plotted according to specified grid parameters.}
\examples{
  # plot generic X, Y and intersection
  plot.new()
  plotVenn2d(rep("",3), radius=c(1.25,1.25), labels=c("Y","X"), Colors=c("yellow","orange","pink"))
  grid.text(expression(paste(X,intersect(Y))),0.5,0.5)
  grid.text("X",0.25,0.5)
  grid.text("Y",0.75,0.5)
  
  # rotated plot
  y <- c(37,29,6)
  labels <- c("A","B")
  plot.new()
  plotVenn2d(y, labels, Colors=rainbow(3), Title = "This is an example of a 2D Venn plot",
    radius=c(0.85, 1.0), rot=45)
}
\keyword{ hplot }
