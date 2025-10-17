\name{plotcircle}
\alias{plotcircle}
\title{
  adds part of a colored circle to a plot
}
\description{
  adds (part of) a colored circle to a plot;
  an arrow can be drawn at a specified position
}

\usage{
plotcircle(r = 1, ...)
}

\arguments{
  \item{r }{radius of circle.
  }
  \item{... }{arguments passed to function \link{plotellipse}.
  }
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{
# symmetrical axes
emptyplot(c(0, 1))
plotcircle(mid = c(0.5, 0.5), r = 0.25, from = 0, to = 3*pi/2,
           arrow = TRUE, arr.pos = 0.5, col = "red")
# symmetrical
plotellipse(mid = c(0.5, 0.5), rx = 0.2, ry = 0.2, 
           arrow = TRUE, arr.pos = 0.5, col = "blue")


#non-symmetrical axes
emptyplot(c(0, 1), c(0, 2), main = "plotcircle", asp = FALSE)
plotcircle(mid = c(0.5, 0.5), r = 0.25, from = 0, to = 3*pi/2,
           arrow = TRUE, arr.pos = 0.5, col = "red")
plotellipse(mid = c(0.5, 0.5), rx = 0.25, ry = 0.25, 
           arrow = TRUE, arr.pos = 0.5, col = "blue")
}

\details{
  \code{plotcircle} calls \code{plotellipse}, making sure that the figure
  drawn effectively looks like a circle. 
  For graphs that have both axes of equal size, the circle will be equal 
  to the ellipse with equal rx and ry. See second example

  see \code{\link{plotellipse}} for details
}

\seealso{
  \code{\link{plotellipse}} to draw ellipses
}

\keyword{aplot}

