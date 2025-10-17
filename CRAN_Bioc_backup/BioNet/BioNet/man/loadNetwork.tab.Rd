\name{loadNetwork.tab}
\alias{loadNetwork.tab}
\title{Load network from tabular format}
\description{
The function loads a network from a tabular format.
}
\usage{
loadNetwork.tab(file, header=TRUE, directed=FALSE, format=c("graphNEL", "igraph"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{File with network to load.}
  \item{header}{Booelan value whether to include header or not.}
  \item{directed}{Booelan value whether the network is to be directed or not.}
  \item{format}{Output format of the network, either \emph{graphNEL} or \emph{igraph}}
}
\author{Marcus Dittrich}
\seealso{\code{\link{loadNetwork.sif}}}
