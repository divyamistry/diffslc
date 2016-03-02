#' @author Divya Mistry
#' @description we create the network, find largest component, save it in RDS
#' 

require(igraph)

#' Read DIP data
intrxns <- read.table(file = "data/minimal_DIP_interaction_data.csv", header = T, sep = ",", na.strings = "-", stringsAsFactors = F)
#' remove self-interaction. These will generate self-loops
intrxns <- intrxns[-which(intrxns$ID.interactor.A == intrxns$ID.interactor.B) ,]

#' create network
mygraph <- graph.edgelist(el = as.matrix(intrxns[,c(1,4)]), directed = F)
#' split into connected components
all.graphs <- decompose.graph(graph = mygraph, min.vertices = 3)
#' There's only one large component, the first one. The other component is 3 nodes, 2 edges.
#' The little two-edge components are difficult to judge for essentiality using proposed method
#' so I'll only look at the largest one.
connGraph <- all.graphs[[1]]
#' don't need these data members now.
rm(all.graphs)

#' functions for DIP to Affy, and DIP to Essential genes mapping
source('src/pre-prep/DipMapper.R')

#' functions in DipMapper.R require two data frames
#'   - YG_S98Annot and DIPtoUniprotMapping
#' Read Yeast annotation
YG_S98Annot <- read.csv(file = "data/YG_S98.na34.annot.csv", header = T, skip=19, na.string = "---", stringsAsFactors = F)
#' Read DIP-to-Uniprot fetched from Uniprot/SwissProt KB
DIPtoUniprotMapping <- read.delim(file = "data/DIP-to-Uniprot-Matched-IDs.tab", header = T, na.string="", stringsAsFactors = F)
colnames(DIPtoUniprotMapping) <- c("UniProtKB", "DIP", "Entry.name", "Protein.names", "Gene.names", 
                                   "Length", "Gene.names.ORF", "Gene.names.ordered.locus", 
                                   "Gene.names.synonym", "Gene.names.primary", "Cross.ref.RefSeq", 
                                   "Cross.ref.EnsembleFungi", "Cross.ref.GeneID", "Cross.ref.CYGD", "Cross.ref.SGD")

#' assign Affy ids using DipMapper. Based on data from UniProt, there are 99 DIP interactors without useful matches
V(connGraph)$Affy.PSs <- sapply(X = V(connGraph)$name, USE.NAMES = T, simplify = T, FUN = function(x) DipMapper(x))

#' assign essentiality
fullDEG <- read.delim(file = "data/degannotation-e.dat", stringsAsFactors = F, na.string = "-")
V(connGraph)$Essential <- sapply(X = V(connGraph)$name, USE.NAMES = T, simplify = T, FUN = function(x) DipEssential(x))

#' Save the processed clean graph. This is the network with DIP data that has been matched to Affy annotation.
#' Expression data will be added later.
saveRDS(object = connGraph, file = "data/Final_hybrid_network.RDS", compress = TRUE)
