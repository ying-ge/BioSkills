\name{BlandAltmanLeh-package}
\alias{BlandAltmanLeh-package}
\alias{BlandAltmanLeh}
\docType{package}
\title{
(Slightly extended) Bland-Altman plots
BlandAltmanLeh
}
\description{
Bland-Altman Plots for assessing agreement between two methods of clinical
measurement and returning associated statistics. Plots are optionally extended
by confidence intervals as described in "J. Martin Bland, Douglas 
G. Altman (1986): Statistical Methods For Assessing Agreement Between Two Methods Of
Clinical Measurement" but not included in the graphics of that publication.
Either base graphics or ggplot2 can be used.
}
\details{
\tabular{ll}{
Package: \tab BlandAltmanLeh\cr
Type: \tab Package\cr
Version: \tab 0.3.1\cr
Date: \tab 2015-12-23\cr
License: \tab GPL\cr
}
 Bland Altman plots are a diagnostic tool for assessing the agreement between
 two methods of measurement or assessing retest reliability from two measurements.
 This package offers plots in base and ggplot2 graphics as well as detailed
 descriptive statistics, thus supporting the construction of individual plots
 based on Bland Altman plots. Bland and Altman describe a way for constructing  
 confidence intervals. This package computes these confidence intervals and 
 includes them into the plots. It also invents the Sunflower-Bland-Altman plot
 for data with ties.
}
\author{
Bernhard Lehnert

Maintainer: Bernhard K. Lehnert <bernhard.lehnert@uni-greifswald.de>
}
\references{
Bland JM, Altman DG, Statistical Methods For Assessing Agreement Between Two 
Methods Of Clinical Measurement, Lancet, 1986; 307-310.

Altman DG, Bland JM, Measurement in medicine: the analysis of method comparison
studies, The Statistician 1983; 32, 307-317.

Vaz S et al., The Case for Using the Repeatability Coeffcient When Calculating
Test-Retest Reliability, PLOS ONE, Sept. 2013, Vol 8, Issue 9.
}
\keyword{ measurement precision, test-retest reliability }
\seealso{ \code{\link{bland.altman.plot},\link{bland.altman.stats}}
}
\examples{
# simple basic Bland Altman plot
a <- rnorm(40,30,10)
b <- 1.01*a + rnorm(40)
bland.altman.plot(a,b, xlab="mean", ylab="difference")

# to get all the data for further analysis
bland.altman.plot(a,b, xlab="mean", ylab="difference", silent=FALSE)

# to include confidence intervals into the plot
bland.altman.plot(a,b, xlab="mean", ylab="difference", conf.int=.95)

# to plot in ggplot2
bland.altman.plot(a,b, graph.sys="ggplot2")

# to mark ties in a Sunflower-Bland-Altman plot
a <- sample(1:5, 40, replace=TRUE)
b <- rep(c(1,2,3,3,5,5,5,5),5)
bland.altman.plot(a, b, sunflower=TRUE)
 
}
