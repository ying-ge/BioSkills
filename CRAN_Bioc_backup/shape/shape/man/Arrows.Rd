\name{Arrows}
\alias{Arrows}
\title{
  adds arrows with improved arrowhead to a plot
}

\description{
  adds one or more arrows  to a plot; arrowhead shape is either curved,
  a triangle, a circle or simple
}

\usage{
Arrows(x0, y0, x1, y1, code = 2, arr.length = 0.4, 
  arr.width = arr.length/2, arr.adj = 0.5, arr.type = "curved",
  segment = TRUE, col = "black", lcol = col, lty = 1, arr.col = lcol, 
  lwd = 1, arr.lwd = lwd, ...)
}

\arguments{
  \item{x0 }{x-coordinates of points *from* which to draw arrows;
    either one value or a vector.
  }
  \item{y0 }{y-coordinates of points *from* which to draw arrows;
    either one value or a vector.
  }
  \item{x1 }{x-coordinates of points *to* which to draw arrows;
    either one value or a vector.
  }
  \item{y1 }{y-coordinates of points *to* which to draw arrows;
    either one value or a vector.
  }
  \item{code }{integer code determining kind of arrows to draw.
  }
  \item{arr.length }{approximate length of arrowhead, in cm;
    either one value or a vector.
  }
  \item{arr.width }{approximate width of arrowhead, in cm;
    either one value or a vector.
  }
  \item{arr.adj }{0,0.5,1 specifying the adjustment of the arrowhead.
  }
  \item{arr.type }{type of arrowhead to draw, one of "none", "simple", "curved",
    "triangle", "circle", "ellipse" or "T".
  }
  \item{segment }{logical specifying whether or not to draw line segments.
  }
  \item{col }{general line color specification; one value or a vector.
  }
  \item{lcol }{line color specifications; either one value or a vector.
    ignored when arr.type = \code{"simple"} or \code{"T"} - use "col"
  }
  \item{lty }{line type specifications; either one value or a vector.
  }
  \item{arr.col }{color of arrowhead; either one value or a vector.
  }
  \item{lwd }{general line width specification. 
  The default value changed to 1 from version 1.4 (was 2) 
  }
  \item{arr.lwd }{line width of arrowhead.
  }
  \item{... }{arguments passed to lines, segments or \link{Arrowhead} function.
  }
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{
xlim <- c(-5 , 5)
ylim <- c(-10, 10)
plot(0, type = "n", xlim = xlim, ylim = ylim, 
     main = "Arrows,  type = 'curved'")
x0 <- runif(100, xlim[1], xlim[2])
y0 <- runif(100, ylim[1], ylim[2])
x1 <- x0+runif(100, -1, 1)
y1 <- y0+runif(100, -1, 1)
Arrows(x0, y0, x1, y1, arr.length = runif(100), code = 2, 
       arr.type = "curved", arr.col = 1:100, lcol = 1:100)
  
plot(0, type = "n", xlim = xlim, ylim = ylim, 
     main = "Arrows,  type = 'circle'")
x0 <- runif(100, xlim[1], xlim[2])
y0 <- runif(100, ylim[1], ylim[2])
x1 <- x0 + runif(100, -1, 1)
y1 <- y0 + runif(100, -1, 1)
Arrows(x0, y0, x1, y1, arr.length = 0.2, code = 3, 
       arr.type = "circle", arr.col = "grey")
  
plot(0, type = "n", xlim = xlim, ylim = ylim, 
     main = "Arrows,  type = 'ellipse'")
Arrows(x0, y0, x1, y1, arr.length = 0.2, arr.width = 0.5, 
       code = 3, arr.type = "ellipse",  arr.col = "grey")

curve(expr = sin(x), 0, 2*pi+0.25, main = "Arrows")
x  <- seq(0, 2*pi, length.out = 10)
xd <- x + 0.025
Arrows(x, sin(x), xd, sin(xd), type = "triangle", 
       arr.length = 0.5, segment = FALSE)
  
xx <- seq(0, 10*pi, length.out = 1000)
plot(sin(xx)*xx, cos(xx)*xx, type = "l", axes = FALSE, 
     xlab = "", ylab = "", main = "Arrows,  type = 'curved'")
x  <- seq(0, 10*pi, length.out = 20)
x1 <- sin(x)*x
y1 <- cos(x)*x
xd <- x+0.01
x2 <- sin(xd)*xd
y2 <- cos(xd)*xd
Arrows(x1, y1, x2, y2, arr.type = "curved", arr.length = 0.4, 
       segment = FALSE, code = 1, arr.adj = 0.5 )

plot(sin(xx)*xx, cos(xx)*xx, type = "l", axes = FALSE, 
     xlab = "", ylab = "", main = "Arrows,  type = 'T'")
Arrows(x1, y1, x2, y2, arr.type = "T", arr.length = 0.4, 
       code = 1, arr.lwd = 2)


# arguments passed to polygon:

xlim <- c(-5 , 5)
ylim <- c(-10, 10)
plot(0, type = "n", xlim = xlim, ylim = ylim, 
     main = "Arrows,  type = 'curved'")
x0 <- runif(100, xlim[1]-1, xlim[2]+0.5)  # exceeds the x-range
y0 <- runif(100, ylim[1], ylim[2])
x1 <- x0+runif(100, -1, 1)
y1 <- y0+runif(100, -1, 1)
Arrows(x0, y0, x1, y1, arr.length = runif(100), code = 2, 
       arr.type = "curved", arr.col = 1:100, lcol = 1:100, xpd = TRUE)


}

\details{
  \code{x0}, \code{y0}, \code{x1}, \code{y1}, \code{arr.length},
  \code{arr.width}, \code{arr.adj}, \code{lcol}, \code{lty} and
  \code{arr.col} can be a vector, of the same length.
   
  For each 'i', an arrow is drawn between the point '(x0[i], y0[i])'
  and the point '(x1[i],y1[i])'.

  \itemize{
    \item If \code{code}=1 an arrowhead is drawn at '(x0[i],y0[i])'
    \item if \code{code}=2 an arrowhead is drawn at '(x1[i],y1[i])'.
    \item If \code{code}=3 an arrowhead is drawn at both ends of the arrow
    \item unless \code{arr.length} = 0, when no head is drawn.
  }
  \itemize{
    \item If \code{arr.adj} = 0.5 then the centre of the arrowhead is at
      the point at which it is drawn.
    \item \code{arr.adj} = 1 causes the tip of the arrowhead to touch
      the point.
    \item \code{arr.adj} = 2 causes the base of the arrowhead to touch
      the point.
  }
  The type of the arrowhead is set with \code{arr.type} which can take
  the values:
  \itemize{
    \item "simple"  : uses comparable \R function \link{arrows}
    \item "triangle": uses filled triangle
    \item "curved"  : draws arrowhead with curved edges
    \item "circle"  : draws circular head
    \item "ellipse" : draws ellepsoid head
    \item "T"       : draws T-shaped (blunt) head
  }
}

\seealso{
  \code{\link{arrows}} the comparable R function

  \code{\link{Arrowhead}}
}

\keyword{aplot}

