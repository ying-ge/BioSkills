\name{drapecol}
\alias{drapecol}
\title{
  draping colors over a persp plot
}

\description{
  generates color(s) that will appear on the surface facets of a "persp" plot.
}
\usage{
drapecol(A, col = femmecol(100), NAcol = "white", lim = NULL)
}

\arguments{
  \item{A }{matrix with input grid.
  }
  \item{col }{color palette.
  }
  \item{NAcol }{color of \code{NA} elements.
  }
  \item{lim }{The limits of the data; if \code{NULL}, the data range will be chosen.
  }
}
\value{
  a vector of character strings giving the colors in hexadecimal format,
  one for each surface facet.
}

\note{
  This function is inspired by a similar function in package \code{fields},
  unfortunately made unavailable in most recent version of fields
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\examples{
persp(volcano, theta = 135, phi = 30, col = drapecol(volcano),
      main = "drapecol")
persp(volcano, theta = 135, phi = 30, col = drapecol(volcano),
      border = NA, main = "drapecol")
}

\seealso{
  \code{\link{persp}}
}

\keyword{aplot}

