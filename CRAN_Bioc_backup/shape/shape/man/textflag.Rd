\name{textflag}
\alias{textflag}
\title{
  adds a filled rounded rectangular box with a text to a plot
}
\description{
  adds a rectangular box with rounded left and right edges to a plot
}

\usage{
textflag(mid, radx, rady, rx = rady, dr = 0.01, 
         col = femmecol(100), lcol = "white", 
         bcol = lcol, lwd = 2, angle = 0, lab = NULL, 
         leftright = TRUE, tcol = NULL, ...)
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
  \item{col }{fill color of the box; the box will be filled from left to right.
  }
  \item{lcol }{line color surrounding box.
  }
  \item{bcol }{line color to remove the ellipse from the rectangular box.
  }
  \item{tcol }{text color.
  }
  \item{lwd }{line width of line surrounding the box.
  }
  \item{angle }{rotation angle, degrees.
  }
  \item{lab }{one label or a vector string of labels to be added in box.
  }
  \item{leftright }{if \code{TRUE} then coloring is from left to right 
    else the coloring is from bottom to top box (for \code{angle} = 0).
  }
  \item{... }{other arguments passed to function \link{text}.
  }
}
\author{Karline Soetaert <karline.soetaert@nioz.nl>}

\examples{
emptyplot()
textflag(mid = c(0.5, 0.5), radx = 0.5, rady = 0.1, 
         lcol = "white",  lab = "hello", cex = 5, font = 2:3)
 
textflag(mid = c(0.5, 0.15), radx = 0.5, rady = 0.1, 
         rx = 0.3,lcol = "black", lab = "hello 2", cex = 4, 
         font = 2, angle = 20, tcol = "darkblue")
  
textflag(mid = c(0.5, 0.85), radx = 0.5, rady = 0.1, rx = 0.03,
         lcol = "white", lab = "hello 3", cex = 4, font = 2, 
         leftright = FALSE)
}
\details{
  \code{radx} and \code{rady} are the horizontal and vertical radiusses of
  the box; \code{rx} is the horizontal radius of the rounded part.

  Here \code{horizontal} and \code{vertical} denote the position BEFORE
  rotation.

  This function is similar to function \link{roundrect}, except that coloring 
  is from left to right.
}

\keyword{aplot}

