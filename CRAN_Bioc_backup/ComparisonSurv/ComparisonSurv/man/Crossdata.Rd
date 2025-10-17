\name{Crossdata}
\docType{data}
\alias{Crossdata}
\title{The Data with Survival Curves Crossed}
\description{
Crossdata is a dataset from simulation, which is violated the proportional hazard ratio assumption.
}
\usage{
data(Crossdata)
}
\format{
A data frame with 200 observations on the following 3 variables.
\describe{
\item{time}{The follow up time for right censored data.}
\item{status}{The status indicator, normally 1=event, 0=alive or right censored.}
\item{group}{The group indicator for comparison, 0= control group, 1= treatment group.}
}
}
\examples{
data(Crossdata)
}
