\name{Descriptive.stat}
\alias{Descriptive.stat}
\title{Summarize Descriptive Statistics Commonly Used in Survival Data}
\description{
A function used to produce result summaries of the results of various descriptive statistics.
}
\usage{
Descriptive.stat(time, status, group, tau="observed", alpha=0.05)
}
\arguments{
\item{time}{The follow up time for right censored data.}
\item{status}{The status indicator, normally 1=event, 0=alive or right censored.}
\item{group}{The group indicator for comparison, and the elements of this vector must take either 0 or 1. Normally, 0= control group, 1= treatment group.}
\item{tau}{The truncation time point (numeric) for the RMST calculation, which needs to be smaller than or equal to the minimum of the largest observed time in each of the two groups. When tau="observed" (default), the minimum of the largest observed time in each of the two groups is used. It can be changed tp tau="event" to choose the minimum of the largest event time in each of the two groups.}
\item{alpha}{The significance level, the default is 0.05.}
}
\value{
A list containing the following components:
\item{result.summary}{The summary statistics in each group, including sample.size, the number of event (event.num), censoring.rate, max.observed.time, max.event.time.}
\item{result.mean}{The results of mean survival time in each group.}
\item{result.quantile}{Quantiles and 0.95 CI of survival time in each group.}
\item{tau}{The truncation time point for the restricted mean survival time (RMST) calculation.}
\item{result.RMST}{The results of RMST in each group. See references for details.}
}
\examples{
#get 'Crossdata' from package
data(Crossdata)
data2<-Crossdata
#
#get result of descriptive statistics with tau=the minimum of the largest observed time
Descriptive.stat(data2$time, data2$status, data2$group,tau="observed")
#get result of descriptive statistics with tau=the minimum of the largest event time
Descriptive.stat(data2$time, data2$status, data2$group,tau="event")
}
