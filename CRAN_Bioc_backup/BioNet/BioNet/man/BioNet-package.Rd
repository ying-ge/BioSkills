\name{BioNet-package}
\Rdversion{1.13.2}
\alias{BioNet-package}
\alias{BioNet}
\docType{package}
\title{
Routines for the functional analysis of biological networks
}
\description{
This package provides functions for the integrated analysis of biological 
networks and the detection of functional modules.
Different datasets can be integrated into the network by assigning p-values derived from statistical 
tests to the nodes of the network. E.g. p-values obtained from the differential expression of
genes from an Affymetrix array are assigned to the nodes of an protein-protein interaction network. 
By fitting a beta-uniform mixture model and calculating scores from the p-values, overall scores of 
network regions can be calculated and an integer linear programming algorithm identifies
the maximum scoring subnetwork.
}
\details{
\tabular{ll}{
Package: \tab BioNet\cr
Type: \tab Package\cr
Version: \tab 1.29.1\cr
Date: \tab 2015-09-11\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}
}
\author{
Marcus Dittrich, Daniela Beisser

Maintainer: Marcus Dittrich <marcus.dittrich@biozentrum.uni-wuerzburg.de>
}
\references{M. T. Dittrich, G. W. Klau, A. Rosenwald, T. Dandekar and T. Mueller (2008) Identifying functional modules in protein-protein interaction networks: an integrated exact approach. \emph{(ISMB2008) Bioinformatics} 24: 13. i223-i231 Jul.


D. Beisser, G. W. Klau, T. Dandekar, T. Mueller and M. Dittrich (2010) BioNet: an R-package for the Functional Analysis of Biological Networks. \emph{Bioinformatics} 26:08. 1129-1130 Apr.
}

