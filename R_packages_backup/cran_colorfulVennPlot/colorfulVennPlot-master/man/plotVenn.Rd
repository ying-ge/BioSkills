\name{plotVenn}
\alias{plotVenn}
\title{Plot Venn diagram for 2-, 3-, 4- dimensional data.}
\description{Plots a 2-, 3-, or 4-dimensional Venn plot depending on the number of input values. The user can specify values, labels for each circle-group and colors. If the vector is of length 3, a 2-d plot is produced. If the vector is of length 7, a 3-d plot is produced. If the vector is length 15, a 4-d plot is produced.
}
\usage{
plotVenn(x, ...)
}
\arguments{
  \item{x}{ for a 3-d plot, a numeric vector of length 7, with a permutation of the
    names c("001","010","011","100","101","110","111").}
  \item{...}{ Additional arguments accepted by plotVenn2d, plotVenn3d or plotVenn4d.}
}
\details{
The '00', '000', '0000' groups are not plotted, so 'percent' data will not add up to 100 percent on the graph.
}
\references{
  Original code: http://elliotnoma.wordpress.com/2011/02/09/venn-diagram/

  Output file: http://elliotnoma.files.wordpress.com/2011/02/venndiagram.png
}
\author{Elliot Noma, Aliona Manvae
}
\note{ It requires package grid to be installed, and can be plotted according to specified grid parameters.}
\examples{
  y <- c(37,29,6,232,121,77,25)
  names(y) <- c("001","010","011","100","101","110","111")
  labels <- c("A","B","C")
  plot.new()
  plotVenn(y, labels, Colors=rainbow(7))
}
\keyword{ hplot }
