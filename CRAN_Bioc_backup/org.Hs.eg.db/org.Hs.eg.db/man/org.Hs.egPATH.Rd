\name{org.Hs.egPATH}
\alias{org.Hs.egPATH}
\alias{org.Hs.egPATH2EG}
\title{Mappings between Entrez Gene identifiers and KEGG pathway identifiers}
\description{
  KEGG (Kyoto Encyclopedia of Genes and Genomes) maintains pathway data
  for various organisms. org.Hs.egPATH maps entrez gene identifiers to the
  identifiers used by KEGG for pathways
}
\details{
  Each KEGG pathway has a name and identifier. Pathway name for a given pathway
  identifier can be obtained using the KEGG data package that can either be
  built using AnnBuilder or downloaded from Bioconductor
  \url{http://www.bioconductor.org}.

  Graphic presentations of pathways are searchable at
  url http://www.genome.ad.jp/kegg/pathway.html by using pathway identifiers as
  keys.   

  Mappings were based on data provided by: KEGG GENOME  
  ftp://ftp.genome.jp/pub/kegg/genomes  
  With a date stamp from the source of: 2011-Mar15

}

\references{
  \url{http://www.genome.ad.jp/kegg/}
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
x <- org.Hs.egPATH
# Get the entrez gene identifiers that are mapped to a KEGG pathway ID
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
  # Get the PATH for the first five genes
  xx[1:5]
  # Get the first one
  xx[[1]]
}

# For the reverse map:
# Convert the object to a list
xx <- as.list(org.Hs.egPATH2EG)
# Remove pathway identifiers that do not map to any entrez gene id
xx <- xx[!is.na(xx)]
if(length(xx) > 0){
    # The entrez gene identifiers for the first two elements of XX
    xx[1:2]
    # Get the first one
    xx[[1]]
}
}
\keyword{datasets}

