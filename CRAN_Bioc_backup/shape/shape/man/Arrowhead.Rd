\name{Arrowhead}
\alias{Arrowhead}
\title{adds arrowheads to a plot}
\description{
  adds one or more arrowheads to a plot; shape is either curved,
  a triangle, a circle or ellipse.
}

\usage{
Arrowhead(x0, y0, angle = 0, arr.length = 0.4, 
          arr.width = arr.length/2, arr.adj = 0.5, 
          arr.type = "curved", lcol = "black", lty = 1, 
          arr.col = lcol, arr.lwd = 2, npoint = 5, ...)
}

\arguments{
  \item{x0 }{x-coordinates of points at which to draw arrowhead;
    either one value or a vector.
  }
  \item{y0 }{y-coordinates of points at which to draw arrowhead;
    either one value or a vector.
  }
  \item{angle }{angle of arrowhead (anti-clockwise, relative to x-axis),
    in degrees [0,360]; either one value or a vector.
  }
  \item{arr.length }{approximate length of arrowhead, in cm;
    either one value or a vector.
  }
  \item{arr.width }{approximate width of arrowhead, in cm;
    either one value or a vector.
  }
  \item{arr.adj }{0,0.5,1 specifying the adjustment of the arrowhead.
  }
  \item{arr.type }{type of arrowhead to draw, one of
    "curved", "triangle", "circle", "ellipse".
  }
  \item{lcol }{line color specifications; either one value or a vector.
  }
  \item{lty }{line type specifications; either one value or a vector.
  }
  \item{arr.col }{color of arrowhead; either one value or a vector.
  }
  \item{arr.lwd }{line width of arrowhead.
  }
  \item{npoint }{only if \code{arr.type} = "curved": number of points to
    draw the curve; increase for smoother arrowheads}
  \item{... }{arguments passed to the \code{polygon} function.
  }}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{
emptyplot(main = "Arrowhead")
Arrowhead(x0 = runif(10), y0 = runif(10), angle = runif(10)*360, 
          arr.length = 0.3, arr.type = "circle", arr.col = "green")
Arrowhead(x0 = runif(10), y0 = runif(10), angle = runif(10)*360, 
          arr.length = 0.4, arr.type = "curved", arr.col = "red")
Arrowhead(x0 = runif(10), y0 = runif(10), angle = runif(10)*360, 
          arr.length = runif(10), arr.type = "triangle", 
          arr.col = rainbow(10))
}

\details{
  \code{x0}, \code{y0}, \code{angle}, \code{arr.length}, \code{arr.width},
  \code{lcol}, \code{lty} and \code{arr.col} can be a vector,
  of the same length.

  \itemize{
     \item if \code{arr.adj} = 0.5, then the centre of the arrowhead is
       at the point at which it is drawn.
     \item \code{arr.adj} = 1 causes the tip of the arrowhead to
       touch the point.
     \item \code{arr.adj} = 0 causes the base of the arrowhead to
       touch the point.
  }
  The type of the arrowhead is set with \code{arr.type} which can take
  the values:
  \itemize{
    \item "triangle": uses filled triangle
    \item "curved"  : draws arrowhead with curved edges
    \item "circle"  : draws circular head (where arr.width=arr.length)
    \item "ellipse" : draws ellipsoid head 
  }
}

\seealso{
  \code{\link{Arrows}}
}

\keyword{aplot}

