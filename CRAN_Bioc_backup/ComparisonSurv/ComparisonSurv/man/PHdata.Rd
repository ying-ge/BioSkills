\name{PHdata}
\alias{PHdata}
\title{The Data Satisfied Proportional Hazard Assumption}
\description{
PHdata is a dataset from simulation, which is satisfied the proportional hazard ratio assumption.
}
\usage{
data(PHdata)
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
data(PHdata)
}



