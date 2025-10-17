#' 16S microbiome data: most diverse bodysites from HMP
#' 
#' The 16S data from the Human Microbiome Project includes only the most
#' diverse bodysites: Antecubital fossa (skin), Stool and Subgingival plaque
#' (oral) and can be analysed using a multilevel approach to account for
#' repeated measurements using our module mixMC. The data include 162 samples
#' (54 unique healthy individuals) measured on 1,674 OTUs.
#' 
#' The data were downloaded from the Human Microbiome Project (HMP,
#' http://hmpdacc.org/HMQCP/all/ for the V1-3 variable region).  The original
#' data contained 43,146 OTU counts for 2,911 samples measured from 18
#' different body sites. We focused on the first visit of each healthy
#' individual and focused on the three most diverse habitats. The prefiltered
#' dataset included 1,674 OTU counts. We strongly recommend to use log ratio
#' transformations on the \code{data.TSS} normalised data, as implemented in
#' the PLS and PCA methods, see details on \url{ www.mixOmics.org/mixMC}.
#' 
#' The \code{data.raw} include a 1 offset in order to be log ratios transformed
#' after TSS normalisation. Consequently, the \code{data.TSS} are TSS
#' normalisation of \code{data.raw}. The CSS normalisation was performed on the
#' orignal data (including zero values)
#' 
#' @name diverse.16S
#' @docType data
#' @usage data(diverse.16S)
#' @format A list containing two data sets, \code{data.TSS} and \code{data.raw}
#' and some meta data information: \describe{ \item{list("data.TSS")}{data
#' frame with 162 rows (samples) and 1674 columns (OTUs). The prefiltered
#' normalised data using Total Sum Scaling normalisation.}
#' \item{list("data.raw")}{data frame with 162 rows (samples) and 1674 columns
#' (OTUs). The prefiltered raw count OTU data which include a 1 offset (i.e. no
#' 0 values).} \item{list("taxonomy")}{data frame with 1674 rows (OTUs) and 6
#' columns indicating the taxonomy of each OTU.} \item{list("indiv")}{data
#' frame with 162 rows indicating sample meta data.}
#' \item{list("bodysite")}{factor of length 162 indicating the bodysite with
#' levels "Antecubital_fossa", "Stool" and "Subgingival_plaque".}
#' \item{list("sample")}{vector of length 162 indicating the unique individual
#' ID, useful for a multilevel approach to taken into account the repeated
#' measured on each individual.} }
#' @return none
#' @references Lê Cao K.-A., Costello ME, Lakis VA, Bartolo, F,Chua XY,
#' Brazeilles R, Rondeau P. MixMC: Multivariate insights into Microbial
#' Communities. PLoS ONE, 11(8): e0160169 (2016).
#' @source The raw data were downloaded from
#' \url{http://hmpdacc.org/HMQCP/all/}. Filtering and normalisation described
#' in our website \url{ www.mixOmics.org/mixMC}
#' @keywords datasets
NULL
