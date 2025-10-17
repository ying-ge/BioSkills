\name{pathifier-package}
\alias{pathifier-package}
\alias{pathifier}
\docType{package}
\title{
Quantify deregulation of pathways in cancer
}
\description{
Pathifier is an algorithm that infers pathway deregulation scores for each tumor sample on the basis of expression data. This score is determined, in a context-specific manner, for every particular dataset and type of cancer that is being investigated. The algorithm transforms gene-level information into pathway-level information, generating a compact and biologically relevant representation of each sample.
}
\details{
\tabular{ll}{
Package: \tab pathifier\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-03-15\cr
License: \tab Artistic-1.0\cr
}
}
\author{
Yotam Drier <drier.yotam@mgh.harvard.edu>
Maintainer: Assif Yitzhaky <assif.yitzhaky@weizmann.ac.il>
}
\references{
Drier Y, Sheffer M, Domany E. Pathway-based personalized analysis of cancer. \emph{Proceedings of the National Academy of Sciences}, 2013, vol. 110(16) pp:6388-6393. (www.pnas.org/cgi/doi/10.1073/pnas.1219651110)

See more information on : http://www.weizmann.ac.il/pathifier/
}

\keyword{ package }

\examples{
data(KEGG) # Two pathways of the KEGG database 
data(Sheffer) # The colorectal data of Sheffer et al.
PDS<-quantify_pathways_deregulation(sheffer$data, sheffer$allgenes,
  kegg$gs, kegg$pathwaynames, sheffer$normals, attempts = 100,
  logfile="sheffer.kegg.log", min_exp=sheffer$minexp, min_std=sheffer$minstd)
}
