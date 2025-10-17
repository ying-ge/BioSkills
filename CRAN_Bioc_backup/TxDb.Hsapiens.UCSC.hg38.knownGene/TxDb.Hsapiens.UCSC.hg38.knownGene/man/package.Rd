\name{TxDb.Hsapiens.UCSC.hg38.knownGene}
\docType{package}

\alias{TxDb.Hsapiens.UCSC.hg38.knownGene-package}
\alias{TxDb.Hsapiens.UCSC.hg38.knownGene}
\alias{TxDb.Hsapiens.UCSC.hg38.knownGene}


\title{Annotation package for TxDb object(s)}

\description{  
  This package loads one or more TxDb objects.  Such TxDb
  objects are an R interface to prefabricated databases contained by
  this package.

  The names of any objects exposed by this package indicate the origin and
  resources exposed.  So for example TxDb.Hsapiens.UCSC.hg19.knownGene
  would be a TxDb object, of Homo sapiens data from UCSC build
  hg19 based on the knownGene Track.
}

\note{
  This data package was made from resources at UCSC on
  2023-09-20 17:25:17 +0000 (Wed, 20 Sep 2023) and based on the hg38 genome based on the knownGene table
}

\author{Bioconductor Core Team, Bioconductor Package Maintainer <maintainer@bioconductor.org> [cre]}


\seealso{
  \link[GenomicFeatures]{transcripts},
  \link[GenomicFeatures]{transcriptsBy}
}

\examples{
## load the library
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
## list the contents that are loaded into memory
ls('package:TxDb.Hsapiens.UCSC.hg38.knownGene')
## show the db object that is loaded by calling it's name
TxDb.Hsapiens.UCSC.hg38.knownGene

}

\keyword{package}
\keyword{data}
