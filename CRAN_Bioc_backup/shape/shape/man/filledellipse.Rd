\name{filledellipse}
\alias{filledellipse}
\title{
  adds a colored and rotated ellipse to a plot
}

\description{
  plots (part of) outer and inner ellipses and colors inbetween;
  color can be a palette
}

\usage{
filledellipse(rx1 = 1, rx2 = 0, ry1 = rx1, ry2 = NULL, mid = c(0,0), 
  dr = 0.01, angle = 0, from = -pi, to = pi, col = femmecol(100),  
  values = NULL, zlim = NULL, lwd = 2, lcol = NA, ...)
}
\value{
 returns, as \code{invisible} a \code{list} containing "xyouter" and "xyinner", the points that define the outer and inner ellipse.
}
\arguments{
  \item{rx1 }{horizontal radius of outer ellipse.
  }
  \item{rx2 }{horizontal radius of inner ellipse.
  }
  \item{ry1 }{vertical radius of outer ellipse.
  }
  \item{ry2 }{vertical radius of inner ellipse.
  }
  \item{mid }{midpoint of ellipse.
  }
  \item{dr }{size of segments, in radians, to draw ellipse
    (decrease for smoother).
  }
  \item{angle }{rotation angle, degrees.
  }
  \item{from }{starting angle for ellipse segment, radians.
  }
  \item{to }{final angle for ellipse segment, radians. The segment is drawn 
    counterclockwise. The default is draw a full ellipse.
  }
  \item{col }{color palette to be used; also allowed are two extremes
    or one value.
  }
  \item{values }{if not \code{NULL}, a matrix providing (radius,z-values)
    couples, used for coloring. .
  }
  \item{zlim }{Only if \code{values} is not \code{NULL}: the minimum and
    maximum z values for which colors should be plotted, defaulting to the
    range of the finite values of the second column of \code{values}.
  }
  \item{lwd }{width of external line.
  }
  \item{lcol }{line color.
  }
  \item{... }{arguments passed to R-function \link{polygon}.
  }
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{
color <- greycol(50)
dr    <- 0.05
emptyplot(xlim = c(-2, 2), ylim = c(-2, 2), col = color[length(color)], 
          main = "filledellipse")
filledellipse(rx1 = 1, mid = c(1, 1) , dr = dr, 
              col = shadepalette(endcol = "darkblue")) 
filledellipse(rx1 = 1, ry1 = 0.5, mid = c(-1, -1), dr = dr, angle = 90, 
              col = shadepalette(endcol = "darkred"))
filledellipse(rx1 = 1, ry1 = 0.5, rx2 = 0.5, dr = dr, mid = c(0, 0), 
              col = c(rev(color), color))
filledellipse(rx1 = 0.5, mid = c(1, -1), dr = dr, from = pi,  to = 1.5*pi, 
              col = rev(shadepalette(endcol = "black"))) 
filledellipse(mid = c(-1, 1))  
  
emptyplot(xlim = c(-2, 2), ylim = c(-2, 2), main = "filledellipse")
filledellipse(rx1 = 0.75, mid = c(-1, 1), col = greycol(100) , dr = dr,  
              values = cbind (1:100, (1:100)^0.5))
filledellipse(rx1 = 0.75, mid = c(1, 1), col = greycol(100)  , dr = dr,  
              values = cbind (1:100, (1:100)))
filledellipse(rx1 = 0.75, mid = c(-1, -1), col = greycol(100), dr = dr,  
              values = cbind (1:100, (1:100)^2))
filledellipse(rx1 = 0.75, mid = c(1, -1), col = greycol(100) , dr = dr,  
              values = cbind (1:100, (1:100)^5))
}

\details{
  draws (part of) an outer and inner ellipse, as specified by inner and
  outer radiusses:

  \code{rx1},\code{ry1}: horizontal and vertical radiusses of outer ellipse;
  \code{rx2},\code{ry2}: same for inner ellipse.
  Here "horizontal" and "vertical" denote the position BEFORE rotation
  
  Fills with a palette of colors inbetween
  
  \code{values}: if not \code{NULL}, a matrix providing (radius,z-values)
  couples, used for coloring.
  Here radius are positive values denoting the relative distance between
  the shapes centre and edge.
  The radiusses are rescaled to be in [0,1] if needed.
  z-values (2nd column of \code{values}) together with \code{zlim} and
  \code{col} denote the coloration level.

  Colors in \code{col} will be interpolated to the z-values and used to
  color an interval as given by the input radiusses.

  If \code{rx2}, the radius of the inner ellipse is 0, the ellipse is full.
}

\seealso{
  \code{\link{filledshape}}, \code{\link{filledcylinder}}
}

\keyword{aplot}

