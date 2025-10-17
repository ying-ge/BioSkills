\name{org.Hs.egORGANISM}
\alias{org.Hs.egORGANISM}
\title{The Organism for org.Hs.eg}
\description{
  org.Hs.egORGANISM is an R object that contains a single item: a
  character string that names the organism for which org.Hs.eg was built.
}

\details{
  Although the package name is suggestive of the organism for which it
  was built, org.Hs.egORGANISM provides a simple way to programmatically
  extract the organism name.
}

\seealso{
  \itemize{
    \item \code{\link[AnnotationDbi]{AnnotationDb-class}} for use of
          the \code{select()} interface.
  }
}

\examples{
## select() interface:
## Objects in this package can be accessed using the select() interface
## from the AnnotationDbi package. See ?select for details.

## Bimap interface:
org.Hs.egORGANISM
}

\keyword{datasets}

