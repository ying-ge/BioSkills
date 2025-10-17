\name{compareNetworks}
\alias{compareNetworks}
\title{
Compare parameters of two networks
}
\description{
The function compares the following parameters of two networks:
diameter, average degree, degree exponent, average path length and plots the
cumulative degree distributions. The networks have to be connected components.
}
\usage{
compareNetworks(network1, network2, plot=TRUE)
}
\arguments{
  \item{network1}{
Network \emph{graphNEL} or \emph{igraph} format.
}
  \item{network2}{
Second network in \emph{graphNEL} or \emph{igraph} format, or subnetwork drawn from first network.
}
  \item{plot}{
Boolean value, whether to plot the cumulative degree distributions.
}
}
\value{
A vector of network parameters is returned:
\item{diam.network1}{Network diameter}
\item{diam.network2}{Diameter of the subnetwork}
\item{av.degree.network1}{Average degree of the network}
\item{av.degree.network2}{Average degree of the subnetwork}
\item{degree.exponent.network1}{Degree exponent of the network}
\item{degree.exponent.network2}{Degree exponent of the subnetwork}
\item{av.path.length.network1}{Average path lenght of the network}
\item{av.path.length.network2}{Average path length of the subnetwork}
}
\author{
Daniela Beisser
}
\examples{
library(DLBCL)
data(interactome)
subnet1 <- largestComp(subNetwork(nodes(interactome)[1:100], interactome))
subnet2 <- largestComp(subNetwork(nodes(interactome)[101:200], interactome))
compareNetworks(network1=subnet1, network2=subnet2)
}