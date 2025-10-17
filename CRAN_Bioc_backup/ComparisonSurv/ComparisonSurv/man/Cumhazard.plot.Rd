\name{Cumhazard.plot}
\alias{Cumhazard.plot}
\title{Plot for Cumulative Hazard Function}
\description{
A function used to plot the curves of cumulative hazard function.
}
\usage{
Cumhazard.plot(time,status,group,col=c(1,4),lwd=c(1,1),lty=c(1,1)
               ,lab.x="",lab.y="",legend=FALSE,local.x=NULL,local.y=NULL
               ,legend.0="",legend.1="")
}
\arguments{
\item{time}{The follow up time for right censored data.}
\item{status}{The status indicator, normally 1=event, 0=alive or right censored.}
\item{group}{The group indicator for comparison, and the elements of this vector must take either 0 or 1. Normally, 0= control group, 1= treatment group.}
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
\seealso{
\code{\link{Survival.plot}},\code{\link{Hazard.plot}}
}
\examples{
#get 'Crossdata' from package
data(Crossdata)
data1<-Crossdata
#
#get cumulative hazard plot
Cumhazard.plot(data1$time, data1$status, data1$group, lab.x='Time', lab.y='Cumulative hazard'
               ,col=c(1,4), lwd=c(3,3), lty=c(1,1), legend=TRUE, local.x=0, local.y=3
               ,legend.0='group0', legend.1='group1')
}
