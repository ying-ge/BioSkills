\name{org.Mm.egACCNUM}
\alias{org.Mm.egACCNUM}
\alias{org.Mm.egACCNUM2EG}
\title{Map Entrez Gene identifiers to GenBank Accession Numbers}  
\description{ 
  org.Mm.egACCNUM is an R object that contains mappings between 
  Entrez Gene identifiers and GenBank accession numbers.    
}
\details{
  This object is a simple mapping of Entrez Gene identifiers
  \url{https://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene} to all possible GenBank 
  accession numbers.

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
x <- org.Mm.egACCNUM
# Get the entrez gene identifiers that are mapped to an ACCNUM
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
  # Get the ACCNUM for the first five genes
  xx[1:5]
  # Get the first one
  xx[[1]]
}
#For the reverse map ACCNUM2EG:
# Convert to a list
xx <- as.list(org.Mm.egACCNUM2EG)
if(length(xx) > 0){
   # Gets the entrez gene identifiers for the first five Entrez Gene IDs
   xx[1:5]
   # Get the first one
   xx[[1]]
}
}

\keyword{datasets}



