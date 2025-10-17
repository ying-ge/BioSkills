\name{Hazard.plot}
\alias{Hazard.plot}
\title{Plot for Kernel-Based Hazard Function}
\description{
A function used to plot the hazard functions. The hazard functions are estimated from right censored data using kernel-based methods.
}
\usage{
Hazard.plot(time,status,group,max.0=NULL,max.1=NULL
            ,col=c(1,4),lwd=c(1,1),lty=c(1,1)
            ,lab.x="",lab.y="",legend=FALSE
            ,local.x=NULL,local.y=NULL,legend.0="",legend.1="")
}
\arguments{
\item{time}{The follow up time for right censored data.}
\item{status}{The status indicator, normally 1=event, 0=alive or right censored.}
\item{group}{The group indicator for comparison, and the elements of this vector must take either 0 or 1. Normally, 0= control group, 1= treatment group.}
\item{max.0}{A numeric value of the maximum horizontal plot coordinate in control group. The default value is max.0=NULL.}
\item{max.1}{A numeric value of the maximum horizontal plot coordinate in treatment group. The default value is max.0=NULL.}
\item{col}{A vector of integers specifying colors for each group. The default value is c(1,4).}
\item{lwd}{A vector of numeric values for line widths. The default value is c(1,1).}
\item{lty}{A vector of integers specifying line types for each curve. The default value is c(1,1).}
\item{lab.x}{Label given to the x-axis.}
\item{lab.y}{Label given to the y-axis.}
\item{legend}{An indicator variable for adding a legend. If it is present this implies legend=FALSE.}
\item{local.x}{Changed the x-axis scale for the legend if legend=TRUE.}
\item{local.y}{Changed the y-axis scale for adding a legend if legend=TRUE.}
\item{legend.0}{Label given to the control group if legend=TRUE.}
\item{legend.1}{Label given to the treatment group if legend=TRUE.}
}
\value{None}
\note{
In Hazard.plot, the kernel estimates are based on the package "muhaz". Here, the algorithm method is the "same bandwidth for all grid points", and the biquadratic is of the boundary kernel function to be used. More details can be found in package "muhaz", the function of muhaz.
}
\seealso{
\code{\link{Survival.plot}},\code{\link{Cumhazard.plot}}
}
\examples{
#get 'Crossdata' from package
data(Crossdata)
data1<-Crossdata
#
#get kernel smooth hazard plot
Hazard.plot(data1$time, data1$status, data1$group, max.0=NULL,max.1=NULL
            , lab.x='Time', lab.y='Smooth hazard ratio'
            , col=c(1,4), lwd=c(3,3), lty=c(1,1), legend=TRUE, local.x=2.7, local.y=1.25
            , legend.0='group0', legend.1='group1')
}
