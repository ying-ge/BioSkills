\name{org.Mm.egCHR}
\alias{org.Mm.egCHR}
\title{Map Entrez Gene IDs to Chromosomes}
\description{
 org.Mm.egCHR is an R object that provides mappings between entrez gene
 identifiers and the chromosome that contains the gene of interest. 
}
\details{
  Each entrez gene identifier maps to a vector of a chromosome. 
  
  Mappings were based on data provided by: Entrez Gene  
  ftp://ftp.ncbi.nlm.nih.gov/gene/DATA  
  With a date stamp from the source of: 2023-Sep11
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
x <- org.Mm.egCHR
# Get the entrez gene identifiers that are mapped to a chromosome
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
  # Get the CHR for the first five genes
  xx[1:5]
  # Get the first one
  xx[[1]]
}
}
\keyword{datasets}

