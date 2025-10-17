\name{resizeCircles}
\alias{resizeCircles}
\title{Match the areas of the 2 circles and their overlap to the input values for X, Y and the intersection of X and Y.}
\description{Calculates the radii of the two circles relative to a fixed center-to-center distance. The areas of the circles and their overlap will be proportionate the three input values of X alone, Y alone and the intersection of X and Y. 
The method used is described in http://elliotnoma.wordpress.com/2013/03/17/spacing-circles-in-a-2-item-venn-diagram/.
}
\usage{
resizeCircles(x_only, y_only, overlap, standardDistance = sqrt(0.5))
}
\arguments{
  \item{x_only}{ a numeric value representing the relative size of X excluding its intersection with Y. }
  \item{y_only}{ a numeric value representing the relaitive size of Y excluding its intersection with X. }
  \item{overlap}{ a numeric value representing the relative size of the intersection of X and Y. }
  \item{standardDistance}{ a numeric value specifying the fixed center-to-center distance to which the radii are matched. The default value is that used by plotVenn2d}
}
\details{
This function is designed to be used in conjunction with the plotVenn2d function. plotVenn2d plots the circles with a fixed center-to-center distance as set by the standardDistance argument. This routine specifies the radii of the circles to control their overlap relative to the fixed center-to-center distance. To increase or decrease the center-to-center distance in plotVenn2d, use the resizePlot parameter. 

The function assumes that the first two arguments are for X-excluding the intersection and Y-excluding the intersection. The total X amount equals this first and third arguments. The total Y amount equals the second and third arguments.
}
\author{Elliot Noma
}
\examples{
  plot.new()
  plotVenn2d(rep("",3), radius=resizeCircles(20,.5, 3), Title=NULL, 
    resizePlot=0.7, labels=c("",""))
}
\keyword{ manip }
