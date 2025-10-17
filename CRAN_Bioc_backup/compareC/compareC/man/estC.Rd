\name{estC}
\alias{estC}
\title{Estimation of a single overall \emph{C} index}
\description{
	This is a function used to estimate a single overall \emph{C} index
}

\usage{
	estC(timeX, statusX, scoreY)
}

\arguments{
	\item{timeX}{The vector of actual survival time \emph{X}, one survival time for each observation}
	\item{statusX}{The matching vector of event indicator for time \emph{X}, 1 if occured and 0 otherwise}
	\item{scoreY}{The vector of the measured biomarker \emph{Y}, one for each of the same observations, based on which the \emph{C} index is to be estimated}
}
\value{
  The estimated \emph{C} index
}

\author{Le Kang, Weijie Chen}

\references{
Harrell FE, Califf RM, Pryor DB, Lee KL, and Rosati RA. (1982) Evaluating the yield of medical tests. The Journal of the American Medical Association, 247(18), 2543--2546

Pencina MJ and D'Agostino RB. (2004) Overall \emph{C} as a measure of discrimination in survival analysis: model specific population value and confidence interval estimation. Statistics in Medicine, 23(13), 2109--2123

Kang L, Chen W, Petrick NA, and Gallas BD. (2014) Comparing two correlated \emph{C} indices with right-censored survival outcome: a one-shot nonparametric approach. Statistics in Medicine, 34(4), 685--703, doi: 10.1002/sim.6370
}


\seealso{
	\code{\link{vardiffC}},\code{\link{compareC}}
}

\examples{
demo(testC)
}
