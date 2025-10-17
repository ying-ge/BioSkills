\name{filledcircle}
\alias{filledcircle}
\title{
  adds colored circle to a plot
}

\description{
  plots (part of) outer and inner circle and colors inbetween;
  color can be a palette.
}

\usage{
filledcircle(r1 = 1, r2 = 0, mid = c(0,0), dr = 0.01, from = -pi, to = pi, 
  col = femmecol(100), values = NULL, zlim = NULL, lwd = 2, lcol = NA, ...)
}

\arguments{
  \item{r1 }{radius of outer circle.
  }
  \item{r2 }{radius of inner circle.
  }
  \item{mid }{midpoint of circle.
  }
  \item{dr }{size of segments, in radians, to draw circle
    (decrease for smoother).
  }
  \item{from }{starting angle for circle segment, radians.
  }
  \item{to }{final angle for circle segment, radians. The segment is drawn 
    counterclockwise. The default is to draw a full circle.
  }
  \item{col }{color palette to be used; also allowed are two extremes or one
    value.
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
color <-graycol(n = 50)
dr    <- 0.05
emptyplot(xlim = c(-2, 2), col = color[length(color)], 
          main = "filledcircle")
filledcircle(r1 = 1, mid = c(1, 1), dr = dr, 
          col = shadepalette(endcol = "darkblue")) 
filledcircle(r1 = 1, mid = c(-1, -1), dr = dr, 
          col = shadepalette(endcol = "darkred"))
filledcircle(r1 = 1, r2 = 0.5, mid = c(0, 0), dr = dr, 
          col = c(rev(color), color))
filledcircle(r1 = 1, mid = c(1, -1), dr = dr, 
          col = intpalette(c("red", "blue", "orange"), 100)) 
filledcircle(mid = c(-1, 1))
  
emptyplot(main = "filledcircle")

for (i in seq(0, 0.45, 0.05)) 
   filledcircle(r1 = i+0.05, r2 = i, 
                mid = c(0.5, 0.5), col = i*20)
}
\value{
 returns, as \code{invisible} a \code{list} containing "xyouter" and "xyinner", the points that define the outer and inner ellipse.
}

\details{
  see \code{\link{filledellipse}} for details
}

\seealso{
  \code{\link{filledshape}}, \code{\link{filledcylinder}},
  \code{\link{filledellipse}}
}

\keyword{aplot}

