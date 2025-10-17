\name{roundrect}
\alias{roundrect}
\title{
  adds a rounded rectangular box to a plot
}
\description{
  adds a rectangular box with rounded left and right edges to a plot
}

\usage{
roundrect(mid, radx, rady, rx = rady, dr = 0.01, 
          col = "white", lcol = "black", lwd = 2, angle = 0, ...)
}
\arguments{
  \item{mid }{midpoint (x,y) of the box.
  }
  \item{radx }{horizontal radius of the box.
  }
  \item{rady }{vertical radius of the box.
  }
  \item{rx }{radius of rounded part.
  }
  \item{dr }{size of segments, in radians, to draw the rounded line 
    (decrease for smoother).
  }
  \item{col }{fill color of the box.
  }
  \item{lcol }{line color surrounding box.
  }
  \item{lwd }{line width of line surrounding the box.
  }
  \item{angle }{rotation angle, degrees.
  }
  \item{... }{arguments passed to function \link{filledshape}.
  }
}
\author{Karline Soetaert <karline.soetaert@nioz.nl>}

\examples{
emptyplot(c(-0.1, 1.1), main = "roundrect")
for (i in 1:10) 
  roundrect(mid = runif(2), col = i, radx = 0.1, rady = 0.05)
for (i in 1:5)
  roundrect(mid = runif(2), col = greycol(20), radx = 0.05,
            rady = 0.05, angle = runif(1)*360)
}
\details{
  \code{radx} and \code{rady} are the horizontal and vertical radiusses of
  the box; \code{rx} is the horizontal radius of the rounded part.

  Here \code{horizontal} and \code{vertical} denote the position BEFORE
  rotation.
}

\keyword{aplot}

