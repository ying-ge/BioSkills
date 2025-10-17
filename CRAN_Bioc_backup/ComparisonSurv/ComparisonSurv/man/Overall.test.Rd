\name{Overall.test}
\alias{Overall.test}
\title{Statistical Inference Methods for Overall Hypothesis Tests}
\description{
A function used to produce the results of various statistical inference methods for overall hypothesis testing, along with the test result of proportional hazards assumption.
}
\usage{
Overall.test(time, status, group, tau=NULL, nperm=500, seed=12345)
}
\arguments{
\item{time}{The follow up time for right censored data.}
\item{status}{The status indicator, normally 1=event, 0=alive or right censored.}
\item{group}{The group indicator for comparison, and the elements of this vector must take either 0 or 1. Normally, 0= control group, 1= treatment group.}
\item{tau}{The truncation time point for the RMST calculation, which needs to be smaller than or equal to the minimum of the largest observed time in each of the two groups. When tau=NULL, the default value (the minimum of the largest event time in each of the two groups) is used.}
\item{nperm}{The times of permutation in area between survival curves(ABS) permutation test, with default nperm=500.}
\item{seed}{The seed number, with default seed=12345.}
}
\value{
A list containing the following components:
\item{method}{Containing 9 statistical inference methods: log-rank, Gehan-Wilcoxon, Tarone-Ware, Weighted Kaplan-Meier (KM), area between survival curves (ABS), ABS permutation, Two-stage, Squared differences, and restricted mean survival time (RMST). See more details of these methods in references.}
\item{statistic}{The statistics of corresponding methods.}
\item{pvalue}{The test P value of corresponding methods.}
}
\note{
The results of Two-stage is based on the "TSHRC" package, more details can be found in package "TSHRC", the function of twostage.

The results of RMST is based on the "survRM2" package, more details can be found in package "survRM2", the function of rmst2.
}
\references{
[1] Li H, Han D, Hou Y, et al. Statistical inference methods for two crossing survival curves: a comparison of methods. PLoS One, 2015, 10: e116774.

[2] Qiu P, Sheng J. A two-stage procedure for comparing hazard rate functions. Journal of the Royal Statistical Society. Series B, Statistical Methodology, 2008, 70: 191-208.

[3] Huang X, Lyu J, Hou Y, et al. A Nonparametric Statistical Method for Two Crossing Survival Curves. Communications in Statistics - Simulation and Computation, 2020. doi: 10.1080/03610918.2020.1753075.
}
\examples{
#get 'Crossdata' from package
data(Crossdata)
data1<-Crossdata
#
#get result of overall hypothesis testing
Overall.test(data1$time,data1$status,data1$group,nperm=10)
}
