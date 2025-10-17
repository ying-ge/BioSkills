\name{Fixpoint.test}
\alias{Fixpoint.test}
\title{Statistical inference methods for testing at a fixed time point}
\description{
A function used to produce the results of various statistical inference methods for testing at a fixed time point.
}
\usage{
Fixpoint.test(time, status, group, t0)
}
\arguments{
\item{time}{The follow up time for right censored data.}
\item{status}{The status indicator, normally 1=event, 0=alive or right censored.}
\item{group}{The group indicator for comparison, and the elements of this vector must take either 0 or 1. Normally, 0= control group, 1= treatment group.}
\item{t0}{The fixed time point for testing.}
}
\value{
A list containing the following components:
\item{est.g0}{The estimation of survival rates at the fixed timepoint for control group.}
\item{est.g1}{The estimation of survival rates at the fixed timepoint for treatment group.}
\item{test}{The results of statistical inference at the fixed timepoint.}
\item{method}{Containing 5 types of statistical inference methods: naive, log, cloglog, arcsin, and logist. See more details in references.}
\item{t0}{The prespecified fixed time point.}
\item{est}{The survival rate at the prespecified timepoint.}
\item{lower.95.CI}{The lower 95 percent pointwise confidence interval for the survival function.}
\item{upper.95.CI}{The upper 95 percent pointwise confidence interval for the survival function.}
\item{statistic}{The statistics of corresponding methods.}
\item{pvalue}{The test P value of corresponding methods.}
}
\references{
[1]Klein JP, Logan B, Harhoff M, et al. Analyzing survival curves at a fixed point in time. Statistics in Medicine, 2007, 26(24):4505-4519.

[2]Anderson JR, Pike LBC. Approximate Confidence Intervals for Probabilities of Survival and Quantiles in Life-Table Analysis. Biometrics, 1982, 38(2):407-416.
}
\examples{
#get 'Crossdata' from package
data(Crossdata)
data1<-Crossdata
#
#if there exist differences at 6 months and 12 months
Fixpoint.test(data1$time,data1$status,data1$group,t0=0.5)
Fixpoint.test(data1$time,data1$status,data1$group,t0=1)
}
