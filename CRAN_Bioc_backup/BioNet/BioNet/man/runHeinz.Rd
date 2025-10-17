\name{runHeinz}
\alias{runHeinz}
\title{
Start HEINZ 
}
\description{
The function starts HEINZ from command line. The HEINZ folder has to 
include the heinz.py python script and the dhea file. CPLEX has to be installed
and accessible from the computer R runs on.  
}
\usage{
runHeinz(heinz.folder="", heinz.e.file, heinz.n.file, N=TRUE, E=FALSE, diff=-1, n=1)
}
\arguments{
  \item{heinz.folder}{
The folder which contains the heinz.py python script and the dhea file.
}
  \item{heinz.e.file}{
The HEINZ edge input file. See \code{\link{writeHeinzEdges}}
}
  \item{heinz.n.file}{
The HEINZ node input file. See \code{\link{writeHeinzNodes}}  
  }
  \item{N}{ 
Boolean value, whether to run HEINZ on nodes.
}
  \item{E}{ 
Boolean value, whether to run HEINZ on edges. HEINZ can run on both with N and E set to TRUE.
}
  \item{diff}{Difference of suboptimal solutions to optimal solution in haming distance in percent. Parameter is set to -1 for optimal solution.}
  \item{n}{Number of optimal and suboptimal solutions, the standard n=1 delivers only the optimal solution.}
}
\details{
This function starts the integer linear programming algorithm to calculate the optimal scoring subnetwork.
The algorithm might be started in the command line when the CPLEX is installed on another machine. To start
it from command line use: heinz.py -e edge.file.txt -n node.file.txt -E False/True -N False/True.
The results can be loaded with \code{\link{readHeinzTree}}, \code{\link{readHeinzGraph}} as a graph object.
}
\references{
M. T. Dittrich, G. W. Klau, A. Rosenwald, T. Dandekar, T. Mueller (2008) Identifying functional modules in protein-protein interaction networks: an integrated exact approach. \emph{(ISMB2008) Bioinformatics}, 24: 13. i223-i231 Jul.
}
\author{
Daniela Beisser
}
\seealso{
\code{\link{writeHeinzEdges}}, \code{\link{writeHeinzNodes}}, \code{\link{readHeinzTree}},  \code{\link{readHeinzGraph}}
}
