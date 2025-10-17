\name{org.Hs.egUCSCKG}
\alias{org.Hs.egUCSCKG}
\title{This mapping has been deprecated and will no longer be available
  after bioconductor 2.6.  See the details section for how you can live
  without it.  For now, it is a map of UCSC "Known Gene" accession numbers with Entrez Gene identifiers}  
\description{ 
  org.Hs.egUCSCKG is an R object that contains mappings between 
  Entrez Gene identifiers and UCSC "Known Gene" accession numbers.    
}

\details{
  This object is a simple mapping of Entrez Gene identifiers
  \url{https://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene} to UCSC
  "Known Gene" Accession Numbers.  Please be advised that UCSC "Known
  Genes" are not more like putative transcripts transcript identifiers
  than they are gene identifiers.  

  Mappings were based on data provided by: UCSC Genome Bioinformatics (Homo sapiens)  
    
  With a date stamp from the source of: 2023-Aug20
  
  This mapping is based on the very latest build available at UCSC for
  this organism as of March 2010.  2.6 is the last release where you
  can expect it to be here.  The GenomicFeatures package contains
  functionality that replaces the need for this mapping.  Please see the
  GenomicFeatures vignette for details on how to use it.
  
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
x <- org.Hs.egUCSCKG
# Get the entrez gene IDs that are mapped to a UCSC KG ID
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
  # Get the UCSC gene IDs for the first five genes
  xx[1:5]
  # Get the first one
  xx[[1]]
}
}
\keyword{datasets}



