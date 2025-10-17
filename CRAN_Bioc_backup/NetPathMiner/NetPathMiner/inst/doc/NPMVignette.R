## ----echo=FALSE, results="hide", warning=FALSE--------------------------------
suppressPackageStartupMessages({library('NetPathMiner')})

## ----Load_package, echo=TRUE, eval=TRUE, results="hide"-----------------------
library(NetPathMiner)

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  graph <- KGML2igraph(filename = file)
#  graph <- SBML2igraph(filename = file)

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  library(rBiopaxParser)
#  biopax = readBiopax(file)
#  graph <- BioPAX2igraph(biopax = biopax)

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  graph <- KGML2igraph(filename = c(file1, file2))

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  graph <- KGML2igraph(filename = ".")

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  # Extract all MIRIAM identifiers from an SBML file.
#  graph <- SBML2igraph(filename = file, miriam = "all")
#  
#  # Extract only miram.go identifiers from a BioPAX file.
#  graph <- BioPAX2igraph(biopax = biopax, miriam = "go")

## ----echo=FALSE, eval=TRUE, results="hide"------------------------------------
file <- file.path(find.package("NetPathMiner"), "extdata", "hsa00860.xml")

## ----echo=TRUE, eval=FALSE, results="hide"------------------------------------
#  graph <- KGML2igraph(filename = file, parse.as = "signaling")
#  
#  graph <- KGML2igraph(filename = file, parse.as = "signaling",
#  	expand.complexes = TRUE)

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
data("ex_sbml")
graph <- ex_sbml
graph

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
head( V(graph) )

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
head( E(graph) )

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
head( V(graph)[ reactions ] )

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
V(graph)[ "reaction_71850" ]$attr

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
getAttrNames(graph)

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
getAttrStatus(graph, pattern = "^miriam.")

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  require("RCurl")
#  # Fetch uniprot annotation
#  graph <- fetchAttribute(graph, organism = "Homo sapiens", target.attr = "miriam.ncbigene" , source.attr = "miriam.uniprot")
#  
#  # Fetch ChEBI annotation.
#  graph <- fetchAttribute(graph, target.attr = "miriam.chebi", source.attr = "miriam.kegg.compound")

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
rgraph <- makeReactionNetwork(graph, simplify=FALSE)
rgraph

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  rgraph <- simplifyReactionNetwork(rgraph)
#  rgraph <- makeReactionNetwork(graph, simplify=TRUE)

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
# Expand complexes of gene network.
ggraph <- expandComplexes(rgraph, v.attr = "miriam.uniprot",
		keep.parent.attr= c("^pathway", "^compartment"))

# Convert reaction network to gene network.
ggraph <- makeGeneNetwork(rgraph)

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
data(ex_microarray)


## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  # Assign weights to edges.
#  if(require("RCurl") && url.exists( NPMdefaults("bridge.web") ))
#  	rgraph <- fetchAttribute(rgraph, organism = "Homo sapiens",
#  						target.attr = "miriam.affy.probeset",
#  						source.attr = "miriam.uniprot")

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  # This requires an internet connection, and RCurl and ALL packages to be present.
#  # Instead, we will actually use a processed ALL data, where features are converted
#  # to miriam.uniprot annotation. (Next chunk)
#  
#  library(ALL)
#  data(ALL)
#  rgraph <- assignEdgeWeights(microarray = exprs(ALL), graph = rgraph,
#  weight.method = "cor", use.attr="miriam.affy.probeset", y=ALL$mol.bio, bootstrap = FALSE)

## ----echo=FALSE, eval=TRUE----------------------------------------------------
# This is what is evaluated.
data(ex_microarray)
rgraph <- assignEdgeWeights(microarray = ex_microarray, graph = rgraph,
weight.method = "cor", use.attr="miriam.uniprot", y=colnames(ex_microarray), bootstrap = FALSE)

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
rgraph$y.labels
head( E(rgraph)$edge.weights )

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
ranked.p <- pathRanker(rgraph, method = "prob.shortest.path",
	K = 25, minPathSize = 6)

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  pathsample <- samplePaths(rgraph, max.path.length = vcount(rgraph),
#  num.samples = 1000, num.warmup = 10)
#  
#  ranked.p <- pathRanker(rgraph, method = "pvalue",
#  sampledpaths = pathsample ,alpha=0.1)

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
# Get paths as edge IDs.
eids <- getPathsAsEIDs(paths = ranked.p, graph = rgraph)

## ----echo=TRUE, eval=TRUE, results="hide"-------------------------------------
# Convert paths to other networks.
eids <- getPathsAsEIDs(paths = ranked.p, graph = ggraph)

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
# Clustering.
ybinpaths <- pathsToBinary(ranked.p)
p.cluster <- pathCluster(ybinpaths, M = 2)

## ----fig=TRUE, pdf=TRUE, echo=TRUE, eval=TRUE---------------------------------
plotClusters(ybinpaths, p.cluster)

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
p.class <- pathClassifier(ybinpaths, target.class = "BCR/ABL", M = 2)

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  plotClassifierROC(p.class)

## ----fig=TRUE, pdf=TRUE, echo=TRUE, eval=TRUE---------------------------------
plotClusters(ybinpaths, p.class)

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
plotNetwork(rgraph, vertex.color="compartment.name")

## ----fig=TRUE, pdf=TRUE, echo=TRUE, eval=FALSE--------------------------------
#  plotPaths(ranked.p, rgraph)
#  
#  # With clusters
#  plotPaths(ranked.p, graph, path.clusters=p.class)

## ----fig=TRUE, pdf=TRUE, echo=TRUE, eval=TRUE---------------------------------
plotAllNetworks(ranked.p, metabolic.net = graph, reaction.net = rgraph,
		path.clusters=p.class, vertex.label = "", vertex.size = 4)

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  layout.c <- clusterVertexByAttr(rgraph, "pathway", cluster.strength = 3)
#  v.color <- colorVertexByAttr(rgraph, "pathway")
#  plotPaths(ranked.p , rgraph, clusters=p.class,
#  	layout = layout.c, vertex.color = v.color)

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  plotCytoscapeGML(graph, file="example.gml", layout = layout.c,
#  				vertex.size = 5, vertex.color = v.color)

## ----echo=TRUE, eval=TRUE, results="hide"-------------------------------------
getGeneSets(graph, use.attr="compartment", gene.attr="miriam.uniprot")

## ----echo=TRUE, eval=TRUE, results="hide"-------------------------------------
getGeneSetNetworks(graph, use.attr="compartment")

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  graphNEL <- toGraphNEL(graph, export.attr="^miriam.")

