\name{getellipse}
\alias{getellipse}
\title{
  x-y coordinates of ellipse
}

\description{
  calculates x-y values for (part of) an ellipse; the ellipse can be rotated
}

\usage{
getellipse(rx = 1, ry = rx, mid = c(0, 0), dr = 0.01, 
           angle = 0, from = -pi, to = pi)
}

\arguments{
  \item{rx }{long radius of ellipse.
  }
  \item{ry }{short radius of ellipse.
  }
  \item{mid }{midpoint of ellipse.
  }
  \item{dr }{size of segments, in radians, to specify ellipse
    (decrease for smoother).
  }
  \item{angle }{rotation angle, degrees.
  }
  \item{from }{starting angle for ellipse segment, radians.
  }
  \item{to }{final angle for ellipse segment, radians. The segment is generated 
    counterclockwise. The default is draw a full ellipse.
  }
}

\value{
  a 2-column matrix with  x-y values of the ellipse
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{
plot(getellipse(1, from = 0, to = pi/2), type = "l", col = "red",
     lwd = 2, main = "getellipse")
lines(getellipse(0.5, 0.25, mid = c(0.5, 0.5)), type = "l", 
      col = "blue", lwd = 2)
lines(getellipse(0.5, 0.25, mid = c(0.5, 0.5), angle = 45), 
      type = "l", col = "green", lwd = 2)

lines(getellipse(0.2, 0.2, mid = c(0.5, 0.5), from = 0, to = pi/2), 
      type = "l", col = "orange", lwd = 2)
lines(getellipse(0.2, 0.2, mid = c(0.5, 0.5), from = pi/2, to = 0), 
      type = "l", col = "black", lwd = 2)
lines(getellipse(0.1, 0.1, mid = c(0.75, 0.5), from = -pi/2, to = pi/2), 
      type = "l", col = "black", lwd = 2)

emptyplot(main = "getellipse")
col <- femmecol(90)
for (i in seq(0, 180, by = 2)) 
      lines(getellipse(0.5, 0.25, mid = c(0.5, 0.5), angle = i),
            type = "l", col = col[(i/2)+1], lwd = 2)
}

\details{
  \code{rx} and \code{ry} are the horizontal and vertical radiusses of the
  ellipses.
  
  points \code{from} and \code{to} are joined counterclockwise. (this has changed
  since version 1.3.4).
}
\seealso{
  \code{\link{plotellipse}}, \code{\link{filledellipse}}
}

\keyword{dplot}

