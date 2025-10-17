\name{mapByVar}
\Rdversion{1.4}
\alias{mapByVar}
\title{
Select probeset by variance and get PPI ID
}
\description{
The function selects for each gene the probeset with the highest variance and gets the
PPI ID for each gene. The PPI identifier is: gene symbol(Entrez ID). Affymetrix identifiers are mapped to the ENTREZ ID.
}
\usage{
mapByVar(exprSet, network=NULL, attr="geneID", ignoreAFFX=TRUE)
}
\arguments{
  \item{exprSet}{
Affymetrix ExpressionSet.
}
  \item{network}{Network that is used to map the Affymetrix identifiers.}
  \item{attr}{The attribute of the network that is used to map the Affymetrix IDs. The IDs are mapped to the unique Entrez gene IDs, which are by default stored in the "geneID" attribute of the network.}
  \item{ignoreAFFX}{Boolean value, whether to ignore or leave AFFX control genes.}
}
\value{
Expression matrix with one gene (PPI ID) per probeset.
}
\author{
Daniela Beisser
}
\examples{
\dontrun{library(ALL);
data(ALL);
mapped.e.set <- mapByVar(ALL);
mapped.e.set[1:10,];}
}