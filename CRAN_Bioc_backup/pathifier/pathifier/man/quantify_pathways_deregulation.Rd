\name{quantify_pathways_deregulation}
\alias{quantify_pathways_deregulation}
\title{
Quantify deregulation of pathways in cancer
}
\description{
Pathifier is an algorithm that infers pathway deregulation scores for each tumor sample on the basis of expression data. This score is determined, in a context-specific manner, for every particular dataset and type of cancer that is being investigated. The algorithm transforms gene-level information into pathway-level information, generating a compact and biologically relevant representation of each sample.
}
\usage{
quantify_pathways_deregulation(data, allgenes, syms, pathwaynames, normals = NULL, 
ranks = NULL, attempts = 100, maximize_stability = TRUE, logfile = "", samplings = NULL,
min_exp = 4, min_std = 0.4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The n x m mRNA expression matrix, where n is the number of genes and m the number of samples.
}
  \item{allgenes}{
A list of n identifiers of genes.
}
  \item{syms}{
A list of p pathways, each pathway is a list of the genes it contains (as appear in "allgenes").
}
  \item{pathwaynames}{
The names of the p pathways.
}
  \item{normals}{
A list of m logicals, true if a normal sample, false if tumor.
}
  \item{ranks}{
External knowledge on the ranking of the m samples, if exists (to use initial guess)
}
  \item{attempts}{
Number of runs to determine stability.
}
  \item{maximize_stability}{
If true, throw away components leading to low stability of sampling noise.
}
  \item{logfile}{
Name of the file the log should be written to (use stdout if empty).
}
  \item{samplings}{
 A matrix specifying the samples that should be chosen in each sampling attempt, chooses a random matrix if samplings is NULL.
}
  \item{min_exp}{
The minimal expression considered as a real signal. Any values below are thresholded to be min_exp.
}
  \item{min_std}{
The minimal allowed standard deviation of each gene. Genes with lower standard deviation are divided by min_std instead of their actual standard deviation. (Recommended: set min_std to be the technical noise).
}
}

\value{
%%  ~Describe the value returned
\item{scores}{The deregulation scores, the main output of pathifier}
\item{genesinpathway}{The genes of each pathway used to devise its dregulation score}
\item{newmeanstd}{Average standart devaition after omitting noisy components}
\item{origmeanstd}{Originial average standart devaition, before omitting noisy components}
\item{pathwaysize}{The number of components used to devise the pathway score}
\item{curves}{The prinicipal curve learned for every pathway}
\item{curves_order}{The order of the points of the prinicipal curve learned for every pathway}
\item{z}{Z-scores of the expression matrix used to learn prinicpal curve}
\item{compin}{The components not omitted due to noise}
\item{xm}{The average expression over all normal samples}
\item{xs}{The standart devation of expression over all normal samples}
\item{center}{The centering used by the PCA}
\item{rot}{The matrix of variable loadings of the PCA}
\item{pctaken}{The number of principal components used}
\item{samplings}{A matrix specifying the samples that should be chosen in each sampling attempt}
\item{sucess}{Pathways for which a deregulation score was sucessfully computed}
\item{logfile}{Name of the file the log was written to}
}

\references{
Drier Y, Sheffer M, Domany E. Pathway-based personalized analysis of cancer. \emph{Proceedings of the National Academy of Sciences}, 2013, vol. 110(16) pp:6388-6393. (www.pnas.org/cgi/doi/10.1073/pnas.1219651110)

See more information on : http://www.weizmann.ac.il/pathifier/
}
\author{
Yotam Drier <drier.yotam@mgh.harvard.edu>
Maintainer: Assif Yitzhaky <assif.yitzhaky@weizmann.ac.il>
}

\examples{
data(KEGG) # Two pathways of the KEGG database 
data(Sheffer) # The colorectal data of Sheffer et al.
PDS<-quantify_pathways_deregulation(sheffer$data, sheffer$allgenes,
  kegg$gs, kegg$pathwaynames, sheffer$normals, attempts = 100,
  logfile="sheffer.kegg.log", min_exp=sheffer$minexp, min_std=sheffer$minstd)
}