\name{filledmultigonal}
\alias{filledmultigonal}

\title{
  adds a colored and rotated multigonal shape to a plot
}

\description{
  draws and colors a rotated shape with equal-sized vertices ;
  color can be a palette.
}

\usage{
filledmultigonal(mid = c(0, 0), rx = 1, ry = rx, nr = 4, 
                 col = femmecol(100), values = NULL, 
                 zlim = NULL, lwd = 2, lcol = NA, angle = 0, ...)
}
\arguments{
  \item{mid }{midpoint of multigonal.
  }
  \item{rx }{horizontal radius.
  }
  \item{ry }{vertical radius.
  }
  \item{nr }{number of sides.
  }
  \item{col }{
    color palette to be used; also allowed are two extremes or one value.
  }
  \item{values }{if not \code{NULL}, a matrix providing (radius,z-values)
    couples, used for coloring.
  }
  \item{zlim }{Only if \code{values} is not \code{NULL}: the minimum and
    maximum z values for which colors should be plotted,
    defaulting to the range of the finite values of the second column of
    \code{values}.
  }
  \item{lwd }{width of external line.
  }
  \item{lcol }{line color.
  }
  \item{angle }{angle of rotation, in degrees.
  }
  \item{... }{arguments passed to R-function \link{polygon}.
  }
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{
emptyplot(c(-1, 1), main = "filledmultigonal")

filledmultigonal(rx = 0.25, ry = 0.125, nr = 3, mid = c(-0.75, 0.75), 
            angle = 45, col = shadepalette("red", "blue", n = 50))
filledmultigonal(rx = 0.125, ry = 0.25, nr = 3, mid = c(-0.25, 0.75), 
            col = shadepalette("red", "yellow", n = 50))
filledmultigonal(rx = 0.25, ry = 0.25 , nr = 3, mid = c(0.25, 0.75),  
            col = c("red", "orange"))
filledmultigonal(rx = 0.25, ry = 0.25 , nr = 3, mid = c(0.75, 0.75), 
            angle = 90,  col = "red")
 
filledmultigonal(rx = 0.25, ry = 0.25, nr = 4, mid = c(-0.75, 0.25), 
            angle = 0, col = shadepalette("red", "blue", n = 50))
filledmultigonal(rx = 0.25, ry = 0.25, nr = 4, mid = c(-0.25, 0.25),  
            angle = 45, col = shadepalette("red", "blue", n = 50))
filledmultigonal(rx = 0.25, ry = 0.125, nr = 4, mid = c(0.25, 0.25), 
            angle = 0, col = shadepalette("red", "blue", n = 50))
filledmultigonal(rx = 0.25, ry = 0.125, nr = 4, mid = c(0.75, 0.25), 
            angle = 45, col = shadepalette("red", "blue", n = 50))
  
filledmultigonal(rx = 0.25, ry = 0.25, nr = 5, mid = c(-0.75, -0.25),  
            angle = 0, col = shadepalette("darkgreen", "lightgreen", n = 50))
filledmultigonal(rx = 0.25, angle = 0, nr = 5, mid = c(-0.25, -0.25), 
            col = rainbow(50))
filledmultigonal(rx = 0.25, angle = 30, nr = 6, mid = c(0.25, -0.25), 
            col = femmecol(50))
filledmultigonal(rx = 0.25, ry = 0.125, angle = 30, nr = 6, mid = c(0.75, -0.25), 
            col = "black")
  
filledmultigonal(rx = 0.25, col = "darkblue", nr = 7, mid = c(-0.75, -0.75))
filledmultigonal(rx = 0.25, col = "darkblue", nr = 9, mid = c(-0.25, -0.75))
filledmultigonal(rx = 0.25, col = "darkblue", nr = 3.7, mid = c(0.25, -0.75))
filledmultigonal(rx = 0.25, col = "darkblue", nr = 4.5, mid = c(0.75, -0.75))
}

\details{
  Coloration proceeds from midpoint to external edge
  
  \code{rx},\code{ry}: horizontal and vertical radiusses of the shape.
  Here "horizontal" and "vertical" denote the position BEFORE rotation

  \code{values}: if not \code{NULL}, a matrix providing (radius,z-values)
  couples, used for coloring.
  Here radius are positive values denoting the relative distance between
  the shapes centre and edge.
  The radiusses are rescaled to be in [0,1] if needed.
  z-values (2nd column of \code{values}) together with \code{zlim} and
  \code{col} denote the coloration level.

  Colors in \code{col} will be interpolated to the z-values and used to
  color an interval as given by the input radiusses.
}
\value{
 returns, as \code{invisible} a \code{list} containing "xyouter" and "xyinner", the points that define the outer and inner ellipse.
}

\seealso{
  \code{\link{filledrectangle}}, \code{\link{filledshape}},
  \code{\link{filledcylinder}}, \code{\link{filledellipse}}
}

\keyword{aplot}

