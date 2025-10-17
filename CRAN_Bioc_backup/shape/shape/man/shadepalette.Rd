\name{shadepalette}
\alias{shadepalette}
\title{
  color palette inbetween two extremes
}
\description{
  Returns color(s) that are a linear interpolation between two colors

  these colors are suitable for shading shapes
}
\usage{
shadepalette(n = 100, endcol = "red", inicol = "white", 
             interval = c(0.0, 1.0))
}
\arguments{
  \item{n }{number of colors.
  }
  \item{endcol }{final color.
  }
  \item{inicol }{initial color.
  }
  \item{interval }{interval *to* where to interpolate.
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
shadepalette(n = 10, "white", "black")
image(matrix(nrow = 1, ncol = 100, data = 1:100),
      col = shadepalette(100, "red", "blue"), main = "shadepalette")
}
\seealso{
  \code{\link{intpalette}}, \code{\link{grey}}, \code{\link{femmecol}}
  \code{\link[grDevices]{colorRamp}} for comparable R functions.
}
\keyword{color}

