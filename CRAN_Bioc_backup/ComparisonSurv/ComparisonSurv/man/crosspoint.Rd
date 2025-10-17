\name{crosspoint}
\alias{crosspoint}
\title{Produce The Cross Time Point When Two Survival Curves Crossed Each Other}
\description{
A function used to produce the cross time points in two survival curves.
}
\usage{
crosspoint(time,status,group)
}
\arguments{
\item{time}{The follow up time for right censored data.}
\item{status}{The status indicator, normally 1=event, 0=alive or right censored.}
\item{group}{The group indicator for comparison, and the elements of this vector must take either 0 or 1. Normally, 0= control group, 1= treatment group.}
}
\value{
A list containing the following components:
\item{crosspoint}{The crossing time points between two survival curves. When there is no crossing between two survival curves, the program will stop with "There is no crossing exists".}
}
\note{
The crosspoint function only used to find the cross time point in survival curves.
}
\examples{
#get 'Crossdata' from package
data(Crossdata)
data1<-Crossdata
#get cross point
#two survival curves crossing at 1.69 years or so
crosspoint(data1$time,data1$status,data1$group)
#
#get 'PHdata' from package
data(PHdata)
data2<-PHdata
#when there exists no crossing
#the program will stop "There is no crossing exists"
# Not run: crosspoint(data2$time,data2$status,data2$group)
}
