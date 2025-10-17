\name{Short.test}
\alias{Short.test}
\title{Statistical Inference Methods for Testing Before a Specified Time Point}
\description{
A function used to produce the results of statistical inference methods for testing the difference before a specified time point.
}
\usage{
Short.test(time,status,group,t0)
}
\arguments{
\item{time}{The follow up time for right censored data.}
\item{status}{The status indicator, normally 1=event, 0=alive or right censored.}
\item{group}{The group indicator for comparison, and the elements of this vector must take either 0 or 1. Normally, 0= control group, 1= treatment group.}
\item{t0}{The cut-point time, which should larger than the minimum of non-censored time and also less than the maximum of all survival time.}
}
\value{
A list containing the following components:
\item{method}{Containing 2 types of statistical inference methods: partial Weighted Kaplan-Meier and partial log-rank.}
\item{t0}{The prespecified cut-point time.}
\item{statistic}{The statistics of corresponding methods.}
\item{pvalue}{The test P value of corresponding methods.}
}
\note{
The partial Weighted Kaplan-Meier method and partial log-rank method are different from the Weighted Kaplan-Meier method and the log-rank method in overall.test. Here, the "partial" means let all observations right censored at the cut-point time t0.
}
\seealso{
\code{\link{crosspoint}},\code{\link{Long.test}}
}

\examples{
#get 'Crossdata' from package
data(Crossdata)
data1<-Crossdata
#
#get cross point
crosspoint(data1$time,data1$status,data1$group)
#two survival curves crossing at 1.69 years
#
#get the result before cross point (t0=1.69)
Short.test(data1$time,data1$status,data1$group,t0=1.69)
}
