\name{compareC}
\alias{compareC}
\title{Testing the difference between two correlated overall \emph{C} indices}
\description{
	This is a function used to test if the difference in two correlated overall \emph{C} indices is statistically significant
}

\usage{
	compareC(timeX, statusX, scoreY, scoreZ)
}

\arguments{
	\item{timeX}{The vector of actual survival time \emph{X}, one survival time for each observation}
	\item{statusX}{The matching vector of event indicator for time \emph{X}, 1 if occured and 0 otherwise}
	\item{scoreY}{The vector of the first measured biomarker or score \emph{Y}, one for each of the same observations}
	\item{scoreZ}{The vector of the second measured biomarker or score \emph{Z}, one for each of the same observations}
}

\value{
 \item{est.c}{The estimated two \emph{C} indices}
 \item{est.diff_c}{The estimated difference of the two \emph{C} indices, i.e.,\eqn{C_{XY}-C_{XZ}}}
 \item{est.vardiff_c}{The estimated variance of the difference of two \emph{C} indices}
 \item{est.varCxy}{The estimated variance of the \emph{C} index for \emph{scoreY}}
 \item{est.varCxz}{The estimated variance of the \emph{C} index for \emph{scoreZ}}
 \item{est.cov}{The estimated covariance between the two \emph{C} indices for \emph{scoreY} and that for \emph{scoreZ}}
 \item{zscore}{Z score (test statistic) for hypothesis testing}
 \item{pval}{P value for the comparison of two \emph{C} indices}
 
}

\note{Under non-random censoring, the returned \emph{p}-value is statistically invalid as the assumption of random censoring is violated.  }


\author{Le Kang, Weijie Chen}

\references{
Harrell FE, Califf RM, Pryor DB, Lee KL, and Rosati RA. (1982) Evaluating the yield of medical tests. The Journal of the American Medical Association, 247(18), 2543--2546

Pencina MJ and D'Agostino RB. (2004) Overall \emph{C} as a measure of discrimination in survival analysis: model specific population value and confidence interval estimation. Statistics in Medicine, 23(13), 2109--2123

Kang L, Chen W, Petrick NA, and Gallas BD. (2014) Comparing two correlated \emph{C} indices with right-censored survival outcome: a one-shot nonparametric approach. Statistics in Medicine, 34(4), 685--703, doi: 10.1002/sim.6370
}


\seealso{
	\code{\link{estC}},\code{\link{vardiffC}}
}

\examples{
demo(testC)
}
