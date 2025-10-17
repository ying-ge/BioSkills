\name{plotellipse}
\alias{plotellipse}
\title{
  adds part of a colored and rotated ellipse to a plot
}
\description{
  adds (part of) a colored, and rotated ellipse to a plot;
  an arrow can be drawn at a specified position.
}
\usage{
plotellipse(rx = 1, ry = 0.2, mid = c(0,0), dr = 0.01, 
            angle = 0, from = -pi, to = pi, type = "l", lwd = 2,
            lcol = "black", col = NULL, arrow = FALSE, 
            arr.length = 0.4, arr.width = arr.length*0.5, 
            arr.type = "curved", arr.pos = 1, arr.code = 2, 
            arr.adj = 0.5, arr.col = "black", ...)
}

\arguments{
  \item{rx }{long radius of ellipse.
  }
  \item{ry }{short radius of ellipse.
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
  \item{to }{final angle for ellipse segment, radians.
  }
  \item{type }{external line or points; "n" if no line.
  }
  \item{lwd }{width of external line.
  }
  \item{lcol }{line color.
  }
  \item{col }{fill color.
  }
  \item{arrow }{drawing arrowhead yes/no.
  }
  \item{arr.length }{length of arrowhead.
  }
  \item{arr.width }{width of arrowhead.
  }
  \item{arr.type }{type of arrow.
  }
  \item{arr.pos }{position of arrow, 0=start,1=end.
  }
  \item{arr.code }{integer code determining kind of arrows to draw.
  }
  \item{arr.adj }{adjustment of arrow.
  }
  \item{arr.col }{color of arrow head.
  }
  \item{... }{arguments passed to R-function \link{lines}.
  }
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{
emptyplot(c(-1, 1), main = "plotellipse")
plotellipse(rx = 0.8, ry = 0.3, angle = 60, col = "blue")
plotellipse(rx = 1.0, ry = 0.6, angle = 0, from = pi, to = 2*pi,
            arrow = TRUE, arr.pos = seq(0.1, 0.5, by = 0.1),
            arr.col = rainbow(5))
plotellipse(rx = 1.0, ry = 0.6, angle = 30, from = pi, to = 1.2*pi,
            col = "red")
plotellipse(rx = 0.1, ry = 0.6, from = 1.5*pi, to = pi,
            lcol = "orange", mid = c(0.2,0.2))
plotellipse(rx = 0.1, ry = 0.6, angle = 30, from = 1.5*pi, to = pi,
            lcol = "orange", mid = c(0.2,0.2))
}

\details{
  \code{rx} and \code{ry} are the horizontal and vertical radiusses of the
  ellipses.

  The ellipse is drawn from the point defined by \code{from} to the point defined
  as \code{to} which are joined anti-clockwise. 

  if \code{arrow} is \code{TRUE}, an arrow is drawn along the path of the
  ellipse.

  \code{arr.length} and \code{arr.width} set the size of the arrow.
    
  The type of the arrowhead is set with \code{arr.type} which can take
  the values:
  \itemize{
    \item "simple"  : uses comparable \R function \link{arrows}.
    \item "triangle": uses filled triangle.
    \item "curved"  : draws arrowhead with curved edges.
    \item "circle"  : draws circular head.
  }
  \code{arr.pos}, a real value between 0 and 1 gives the position
  (0=start,1=end).

  \code{arr.col} specifies the color, \code{arr.code} specifies where the
  angle points to.

  \code{arr.adj} specifies the position adjustment - see \code{\link{Arrows}}
  for details.
}

\seealso{
  \code{\link{getellipse}}, \code{\link{filledellipse}},
  \code{\link{plotcircle}}.
}

\keyword{aplot}

