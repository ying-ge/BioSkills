\name{femmecol}
\alias{femmecol}
\title{red-green-blue color palette}
\description{
  Creates a vector of (n) contiguous colors
  (darkblue-blue-cyan-yellow-red-darkred).
}

\usage{
femmecol(n = 100)
}

\arguments{
  \item{n }{number of colors.
  }
}

\value{
  a vector of character strings giving the colors in hexadecimal format
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{
filled.contour(volcano, color = femmecol, asp = 1, main = "femmecol")
femmecol(10)
image(matrix(nrow = 1, ncol = 100, data = 1:100), 
      col = femmecol(100), main = "femmecol")
}

\seealso{
  \code{\link{rainbow}}, \code{\link{heat.colors}},
  \code{\link{topo.colors}}, the comparable R-functions.

  \code{\link{intpalette}}, \code{\link{shadepalette}}
}

\keyword{color}

