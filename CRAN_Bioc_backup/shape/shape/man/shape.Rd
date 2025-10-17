\name{shape-package}
\alias{shape-package}
\alias{shape}
\docType{package}
\title{
  Functions for plotting graphical shapes, colors
}
\description{
  Functions for plotting graphical shapes
  such as ellipses, circles, cylinders, arrows, ...
     
  Support for the book "A practical guide to ecological modelling -
  using R as a simulation platform"
  by Karline Soetaert and Peter M.J. Herman (2009). Springer.
}

\details{
  \tabular{ll}{
    Package: \tab shape\cr
    Type: \tab Package\cr
    Version: \tab 1.3.4\cr
    Date: \tab 2011-07-30\cr
    License: \tab GNU Public License 3 or above\cr
  }
This package is used in R-package ecolMod, which includes many more examples.

See also R-package diagram.

Changes in version 1.3.4: more consistent drawing of ellipse and circle segments,
(functions getellipse, getcircle), added textflag. (both suggested by Tom Wilson)
}
\author{
  Karline Soetaert (Maintainer)
}
\seealso{
  \code{\link{A4}}, \code{\link{writelabel}}, \code{\link{emptyplot}},
  \code{\link{drapecol}}, \code{\link{femmecol}}, \code{\link{intpalette}},
  \code{\link{shadepalette}}, \code{\link{colorlegend}},
  \code{\link{greycol}}, \code{\link{rotatexy}}, \code{\link{Arrowhead}},
  \code{\link{Arrows}}, \code{\link{cylindersegment}},
  \code{\link{filledcylinder}}, \code{\link{filledcircle}},
  \code{\link{filledellipse}}, \code{\link{filledmultigonal}},
  \code{\link{filledrectangle}}, \code{\link{filledshape}},
  \code{\link{getellipse}}, \code{\link{plotcircle}},
  \code{\link{plotellipse}}, \code{\link{roundrect}}, \code{\link{textflag}}.
}

\examples{
\dontrun{
## show examples (see respective help pages for details)
example(rotatexy)
example(filledshape)

## run demos
demo("colorshapes") # creating colored shapes

## open the directory with source code of demos
browseURL(paste(system.file(package="shape"), "/demo", sep=""))

## show package vignette
vignette("shape")
edit(vignette("shape"))
browseURL(paste(system.file(package="shape"), "/doc", sep=""))
}
}

\keyword{ package }