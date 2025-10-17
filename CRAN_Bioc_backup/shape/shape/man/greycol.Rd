\name{greycol}
\alias{greycol}
\alias{graycol}
\title{
  white-black color palette
}

\description{
  Creates a vector of (n) contiguous colors from white/grey to black
}

\usage{
greycol(n = 100, interval = c(0.0, 0.7))
}

\arguments{
  \item{n }{number of colors.
  }
  \item{interval }{interval *to* where to interpolate.
  }
}
\value{
  a vector of character strings giving the colors in hexadecimal format.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{
filled.contour(volcano, color = graycol, asp = 1, main = "greycol,graycol")
graycol(10)
image(matrix(nrow = 1, ncol = 100, data = 1:100),
      col = graycol(100), main = "greycol,graycol")
}

\details{
  greycol is an alias of graycol
}

\seealso{
  \code{\link{rainbow}}, \code{\link{heat.colors}},
  \code{\link{topo.colors}}, \code{\link{femmecol}}
}

\keyword{color}

