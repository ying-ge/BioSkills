\name{save3dModule}
\Rdversion{1.4}
\alias{save3dModule}
\title{
Save a 3D plot of the network
}
\description{
The function saves a 3D plot of a network to file, therefore it requires the plot to be open. A screenshot of the 3D plot can be saved in "pdf" format.
Background of the device is changed to white for plotting. The screenshot can take several seconds for large plots.
}
\usage{
save3dModule(file)
}
\arguments{
\item{file}{File to save to.}
}
\author{
Daniela Beisser
}
\seealso{\code{\link{plot3dModule}}, \code{\link{plotModule}}}
\examples{
library(DLBCL)
data(dataLym)
data(interactome)
interactome <- subNetwork(dataLym$label, interactome)
fchange <- dataLym$diff
names(fchange) <- dataLym$label
subnet <- largestComp(subNetwork(nodes(interactome)[1:100], interactome))
diff <- fchange[nodes(subnet)]

\dontrun{library(rgl);
plot3dModule(network=subnet, diff.or.score=diff);
save3dModule(file="test")}
}