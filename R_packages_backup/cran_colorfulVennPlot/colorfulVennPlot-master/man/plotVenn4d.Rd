\name{plotVenn4d}
\alias{plotVenn4d}
\title{Plot Venn diagram for 4-dimensional data.}
\description{Given a vector of 15 values, which describe 4-dimensional data, it plots a Venn diagram, i.e. 'crossing circles'. The user is able to specify values, labels for each circle-group and colors.
}
\usage{
plotVenn4d(x, labels = c('A','B','C','D'),
  Colors = c('red', 'yellow', 'green', 'pink','darkgreen','blue','lightblue','tan', 
  'yellowgreen','orange','purple','white','grey','plum','brown'),
  Title = NULL, shrink = 1, rot=45)
}
\arguments{
  \item{x}{ a numeric vector of length 15, with names c('1000', '0100', '1100', '0010', '1010', '0110', '1110', 
  '0001', '1001', '0101', '1101', '0011', '1011', '0111', '1111'))in suitable order.}
  \item{labels}{ a character vector of length 4, providing names for the 4 dimensions. }
  \item{shrink}{ a numeric value specifying zooming effect of the plot. }
  \item{Colors}{ a vector of color names for the backgrounds of each part of the diagram. }
  \item{Title}{ a string containing the graph title. }
  \item{rot}{ a numeric value for the number of degrees to rotate the graph. }
}
\details{
NOTICE: This only works for 4-dimensional data.

By definition, value '0000' lays outside the plotted diagram. Therefore, it will not be plotted. Because of that, if plotting 'percent' data, all values will not add up to 100 percent, missing 0000's allocation.
}
\author{ Elliot Noma, Aliona Manvae
}
\note{ It requires package grid to be installed, and can be plotted according to specified grid parameters.}
\examples{

plot.new()
Colors <- c('red', 'yellow', 'green', 'pink', 'darkgreen','blue','lightblue','tan', 
  'yellowgreen','orange','purple','white','grey','plum','brown')
regions <- seq(15)
names(regions) <- c('1000', '0100', '1100', '0010', '1010', '0110', '1110', '0001', 
  '1001', '0101', '1101', '0011', '1011', '0111', '1111')
plotVenn4d(regions, Colors=Colors, Title = '4-dimensional')

}
\keyword{ hplot }
