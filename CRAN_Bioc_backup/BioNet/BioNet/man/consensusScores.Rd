\name{consensusScores}
\Rdversion{1.4}
\alias{consensusScores}
\title{
Calculation of a consensus score for a network
}
\description{
The function calculates consensus scores for a network, given a list of replicate modules.
}
\usage{
consensusScores(modules, network, ro=length(modules)/2)
}
\arguments{
  \item{modules}{Calculated modules from pseudo-replicates of expression values in \emph{igraph} or \emph{graphNEL} format.}
  \item{network}{Interaction network, which shoupld be scores. In \emph{igraph} or \emph{graphNEL} format}
  \item{ro}{Threshold which is subtracted from the scores to obtain positive and negative value. The default value is half of the number of replicates.}
}
\value{
A result list is returned, consisting of:
\item{N.scores}{Numerical vector node scores.}
\item{E.scores}{Numerical vector edge scores.}
\item{N.frequencies}{Numerical vector node frequencies from the replicate modules.}
\item{E.frequencies}{Numerical vector edge frequencies from the replicate modules.}
}
\author{
Daniela Beisser
}
\examples{
library(DLBCL)
data(interactome)
network <- interactome
# precomputed Heinz modules from pseudo-replicates
\dontrun{lib <- file.path(.path.package("BioNet"), "extdata") 
modules <- readHeinzGraph(node.file=file.path(datadir, "ALL_n_resample.txt.0.hnz"), network=network)
cons.scores <- consensusScores(modules, network)
}
}
