\name{intpalette}
\alias{intpalette}
\title{color palettes}
\description{
  Returns color(s) that are a linear interpolation of a given set of colors.
}

\usage{
intpalette(inputcol, numcol = length(x.to), x.from = NULL, x.to = NULL)
}

\arguments{
  \item{inputcol }{initial colors, *from* where to interpolate.
  }
  \item{numcol }{number of colors to interpolate *to*.
  }
  \item{x.from }{x-values *from* where to interpolate.
  }
  \item{x.to }{x-values where to interpolate *to*.
  }
}

\value{
  a vector of character strings giving the interpolated colors in
  hexadecimal format
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{
intpalette(c("white", "black"), n = 10)
grey(seq(1, 0, length.out = 10))
image(matrix(nrow = 1,ncol = 100, data = 1:100),
      col = intpalette(c("red", "blue"), numcol = 100),
      main = "intpalette")
image(matrix(nrow = 1, ncol = 100, data = 1:100),
      col = intpalette(c("red", "blue", "yellow"), numcol = 100),
      main = "intpalette")
}

\details{
  Return value is a vector of *colors* in hexadecimal format.
     
  This is different from \code{\link[grDevices]{colorRamp}}(R function),
  that returns a *function*
}

\seealso{
  \code{\link{greycol}}, \code{\link{femmecol}}, \code{\link{shadepalette}},
  \code{\link[grDevices]{colorRamp}} for comparable R function
}

\keyword{color}

