\name{Survival.plot}
\alias{Survival.plot}
\title{Plot for Cumulative Survival Function}
\description{
A function used to plot the survival curves.
}
\usage{
Survival.plot(time,status,group,...)
}
\arguments{
\item{time}{The follow up time for right censored data.}
\item{status}{The status indicator, normally 1=event, 0=alive or right censored.}
\item{group}{The group indicator for comparison, and the elements of this vector must take either 0 or 1. Normally, 0= control group, 1= treatment group.}
\item{...}{Parameters in function plot()}
}
\value{None}
\note{
Survival.plot is based on the survival package, more details can be found in package "survival", the function of plot.survfit.
}
\seealso{
\code{\link{Cumhazard.plot}},\code{\link{Hazard.plot}}
}
\examples{
#get 'Crossdata' from package
data(Crossdata)
data1<-Crossdata
#
#get survival plot
Survival.plot(data1$time, data1$status, data1$group, xlab = "Time", ylab = "Survival"
              , col = c(1, 4), lwd = c(3, 3), lty= c(1, 1))
legend(2.5, 1, c("group0", "group1"), col = c(1, 4), lwd = c(3, 3), lty= c(1, 1), bty="n")
}
