\name{org.Mm.egMGI}
\alias{org.Mm.egMGI}
\alias{org.Mm.egMGI2EG}
\title{Map MGI gene Acession numbers with Entrez Gene identifiers}  
\description{ 
  org.Mm.egMGI is an R object that contains mappings between 
  Entrez Gene identifiers and Jackson Laboratory MGI gene accession numbers.    
}
\details{
  This object is a simple mapping of Entrez Gene identifiers
  \url{https://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene} to MGI
  gene Accession Numbers.

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
x <- org.Mm.egMGI
# Get the entrez gene IDs that are mapped to an MGI ID
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
  # Get the MGI gene IDs for the first five genes
  xx[1:5]
  # Get the first one
  xx[[1]]
}
#For the reverse map MGI2EG:
# Convert to a list
xx <- as.list(org.Mm.egMGI2EG)
if(length(xx) > 0){
   # Gets the entrez gene IDs for the first five MGI IDs
   xx[1:5]
   # Get the first one
   xx[[1]]
}
}
\keyword{datasets}



