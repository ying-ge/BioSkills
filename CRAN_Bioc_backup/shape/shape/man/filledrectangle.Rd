\name{filledrectangle}
\alias{filledrectangle}

\title{
  adds a colored and rotated rectangle to a plot
}

\description{
  plots and colors a rotated rectangle; color can be a palette
}

\usage{
filledrectangle(mid = c(0, 0), wx = 1, wy = wx, col = femmecol(100), 
                values = NULL, zlim = NULL, lwd = 2, lcol = NA, 
                angle = 0, ...)
}

\arguments{
  \item{mid }{midpoint of rectangle.
  }
  \item{wx }{horizontal width.
  }
  \item{wy }{vertical width.
  }
  \item{col }{color palette to be used; also allowed are two extremes or
    one value.
  }
  \item{values }{ if not \code{NULL}, a matrix providing (radius,z-values)
    couples, used for coloring.
  }
  \item{zlim }{ Only if values is not NULL: the minimum and maximum z values
    for which colors should be plotted, defaulting to the range of the
    finite values of the second column of \code{values}.
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
color <- shadepalette(grey(0.3), "lightblue", n = 50)
emptyplot(main = "filledrectangle")
filledrectangle(wx = 0.5, wy = 0.5, col = color, 
                mid = c(0.5, 0.5), angle = 0)
filledrectangle(wx = 0.25, wy = 0.25, col = "darkblue", 
                mid = c(0.5, 0.5), angle = 45)
filledrectangle(wx = 0.125, wy = 0.125, col = c("lightblue","blue"),
                mid = c(0.5, 0.5), angle = 90)
  
color <- shadepalette(grey(0.3), "blue", n = 50)
emptyplot(c(-1, 1), main = "filledrectangle")
filledrectangle(wx = 0.5, wy = 0.5, col = color, 
                mid = c(0, 0), angle = 0)
filledrectangle(wx = 0.5, wy = 0.5, col = color, 
                mid = c(0.5, 0.5), angle = 90)
filledrectangle(wx = 0.5, wy = 0.5, col = color, 
                mid = c(-0.5, -0.5), angle = -90)
filledrectangle(wx = 0.5, wy = 0.5, col = color, 
                mid = c(0.5, -0.5), angle = 180)
filledrectangle(wx = 0.5, wy = 0.5, col = color, 
                mid = c(-0.5, 0.5), angle = 270)
}

\details{
  If \code{angle}=0, coloration starts from top to bottom. This is different
  from \code{\link{filledmultigonal}},
  where coloration proceeds from middle to external
    
  \code{wx},\code{wy}: horizontal and vertical width of the shape
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
\seealso{
  \code{\link{filledmultigonal}}, \code{\link{filledshape}},
  \code{\link{filledcylinder}}, \code{\link{filledellipse}}

  \code{\link{polygon}}, \code{\link{rect}} for corresponding R-functions.
}

\keyword{aplot}

