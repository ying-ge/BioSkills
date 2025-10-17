\name{KEGG}
\alias{KEGG}
\docType{data}
\title{
Two pathways of the KEGG database
}
\description{
Two pathways (MISMATCH REPAIR and REGULATION OF AUTOPHAGY) of the KEGG database
}
\usage{data(KEGG)}
\format{
  \describe{
    \item{\code{pathwaynames}}{The names of the pathways}
    \item{\code{gs}}{The list of genes (by official gene symbol) in each pathway}
  }
}

\source{
Kanehisa M, Goto S, Sato Y, Furumichi M and Tanabe M. KEGG for integration and interpretation of large-scale molecular datasets. \emph{Nucleic Acids Res}, 2012, Vol 40(Database issue):D109-D114.
}
\examples{
data(KEGG)
}
\keyword{datasets}
