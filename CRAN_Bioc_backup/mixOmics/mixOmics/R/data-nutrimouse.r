#' Nutrimouse Dataset
#' 
#' The \code{nutrimouse} dataset contains the expression measure of 120 genes
#' potentially involved in nutritional problems and the concentrations of 21
#' hepatic fatty acids for forty mice.
#' 
#' The data sets come from a nutrigenomic study in the mouse (Martin \emph{et
#' al.}, 2007) in which the effects of five regimens with contrasted fatty acid
#' compositions on liver lipids and hepatic gene expression in mice were
#' considered. Two sets of variables were acquired on forty mice:
#' 
#' \itemize{ \item gene: expressions of 120 genes measured in liver cells,
#' selected (among about 30,000) as potentially relevant in the context of the
#' nutrition study. These expressions come from a nylon macroarray with
#' radioactive labelling;
#' 
#' \item lipid: concentrations (in percentages) of 21 hepatic fatty acids
#' measured by gas chromatography. }
#' 
#' Biological units (mice) were cross-classified according to two factors
#' experimental design (4 replicates):
#' 
#' \itemize{ \item Genotype: 2-levels factor, wild-type (WT) and
#' PPAR\eqn{\alpha} -/- (PPAR).
#' 
#' \item Diet: 5-levels factor. Oils used for experimental diets preparation
#' were corn and colza oils (50/50) for a reference diet (REF), hydrogenated
#' coconut oil for a saturated fatty acid diet (COC), sunflower oil for an
#' Omega6 fatty acid-rich diet (SUN), linseed oil for an Omega3-rich diet (LIN)
#' and corn/colza/enriched fish oils for the FISH diet (43/43/14). }
#' 
#' @name nutrimouse
#' @docType data
#' @usage data(nutrimouse)
#' @format A list containing the following components: \describe{
#' \item{list("gene")}{data frame with 40 observations on 120 numerical
#' variables.} \item{list("lipid")}{data frame with 40 observations on 21
#' numerical variables.} \item{list("diet")}{factor of 5 levels containing 40
#' labels for the diet factor.} \item{list("genotype")}{factor of 2 levels
#' containing 40 labels for the diet factor.} }
#' @return none
#' @references Martin, P. G. P., Guillou, H., Lasserre, F., Déjean, S., Lan,
#' A., Pascussi, J.-M., San Cristobal, M., Legrand, P., Besse, P. and Pineau,
#' T. (2007). Novel aspects of PPAR\eqn{\alpha}-mediated regulation of lipid
#' and xenobiotic metabolism revealed through a multrigenomic study.
#' \emph{Hepatology} \bold{54}, 767-777.
#' @source The \code{nutrimouse} dataset was provided by Pascal Martin from the
#' Toxicology and Pharmacology Laboratory, National Institute for Agronomic
#' Research, French.
#' @keywords datasets
NULL
