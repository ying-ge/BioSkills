\name{Sheffer}
\alias{Sheffer}
\docType{data}
\title{
Sheffer et al. colorectal dataset
}
\description{
Partial data from Sheffer et al. paper
}
\usage{data(Sheffer)}
\format{
  \describe{
    \item{\code{data}}{the expression data}
    \item{\code{samples}}{sample names}
    \item{\code{normals}}{which of the samples is a normal sample}
    \item{\code{minstd}}{minimal standart deviation allowed}
    \item{\code{minexp}}{minimal value of experssion allowed}
    \item{\code{allgenes}}{the list of genes (by official gene symbol)}
  }
}

\source{
Sheffer et.\ al. Association of survival and disease progression with chromosomal instability: A genomic exploration of colorectal cancer.
\emph{PNAS}, 2009, Vol 106(17) pp: 7131-7136.
}
\examples{
data(Sheffer)
}
\keyword{datasets}
