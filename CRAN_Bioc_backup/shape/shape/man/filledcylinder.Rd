\name{filledcylinder}
\alias{filledcylinder}
\title{
  adds a colored and rotated cylinder to a plot
}

\description{
  adds a rotated and colored cylinder to a plot; color can be a palette
}

\usage{
filledcylinder(rx = 1, ry = rx, len = 1, col = femmecol(100), 
              lcol = NA, lwd = 2, lcolint = NULL, ltyint = 1, 
              lwdint = lwd, mid = c(0,0), angle = 0, delt = 1,
              dr = 0.01, topcol = NULL, botcol = NULL, ...)
}

\arguments{
  \item{rx }{horizontal radius.
  }
  \item{ry }{vertical radius.
  }
  \item{len }{length.
  }
  \item{col }{color palette to be used; also allowed are two extremes
    or one value.
  }
  \item{lcol }{line color on external surface.
  }
  \item{lwd }{only if lcol!=NA, width of external line.
  }
  \item{lcolint }{only if lcol!=NA, line color on internal (hidden) surface.
  }
  \item{ltyint }{only if lcol!=NA, line type on internal (hidden) surface.
  }
  \item{lwdint }{only if dlcol!=NA, line width on internal (hidden) surface.
  }
  \item{mid }{midpoint of cylinder.
  }
  \item{angle }{rotation angle, degrees.
  }
  \item{delt }{increase factor, from left to right.
  }
  \item{dr }{size of segments, in radians, to draw top/bottom ellipse 
          (decrease for smoother).
  }
  \item{topcol }{color (palette) of top (right) surface.
  }
  \item{botcol }{color (palette) of bottom (left) surface.
  }
  \item{... }{arguments passed to function \link{filledellipse}.
  }
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{
emptyplot(c(-1.2, 1.2), c(-1, 1), main = "filledcylinder")
col  <- c(rev(greycol(n = 50)), greycol(n = 50))
col2 <- shadepalette("red", "blue", n = 50)
col3 <- shadepalette("yellow", "black", n = 50)
filledcylinder(rx = 0., ry = 0.2, len = 0.25, angle = 0, col = col,  
          mid = c(-1, 0), topcol = col[25])          
filledcylinder(rx = 0., ry = 0.2, angle = 90, col = col, 
          mid = c(-0.5, 0), topcol = col[25])          
filledcylinder(rx = 0.1, ry = 0.2, angle = 90, col = c(col2, rev(col2)), 
          mid = c(0.45, 0), topcol = col2[25])          
filledcylinder(rx = 0.05, ry = 0.2, angle = 90, col = c(col3, rev(col3)), 
          mid = c(0.9, 0), topcol = col3[25])          
filledcylinder(rx = 0.1, ry = 0.2, angle = 90, col = "white",  
          lcol = "black", lcolint = "grey")          
  
emptyplot(c(-1, 1), c(-1, 1), main = "filledcylinder")
col  <- shadepalette("blue", "black", n = 50)
col2 <- shadepalette("red", "black", n = 50)
col3 <- shadepalette("yellow", "black", n = 50)
filledcylinder(rx = 0.025, ry = 0.2, angle = 90, col = c(col2, rev(col2)), 
          mid = c(-0.8, 0), topcol = col2[25], delt = -1, lcol = "black")
filledcylinder(rx = 0.1, ry = 0.2, angle = 00, col = c(col, rev(col)), 
          mid = c(0.0, 0.0), topcol = col, delt = -1.2, lcol = "black")
filledcylinder(rx = 0.075, ry = 0.2, angle = 90, col = c(col3, rev(col3)), 
          mid = c(0.8, 0), topcol = col3[25], delt = 0.0, lcol = "black")
}

\details{
  When angle = 0 (the default),  the cylinder is parallel to the x-axis
  
  \code{rx} and \code{ry} are the horizontal and vertical radiusses of the
  bordering ellipses.
  Here "horizontal" and "vertical" denote the position BEFORE rotation

  if \code{delt} > 1, the width of the cylinder will increase from left
  to right.
}

\seealso{
  \code{\link{filledellipse}}, \code{\link{filledshape}}
}

\keyword{aplot}

