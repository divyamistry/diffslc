#' @author Divya Mistry
#' @description Network statistics for hybrid PIN with GSE4341
#'
if (!is.element('pROC', installed.packages()[,1])) { stop("Install pROC package.") }
if (!is.element('igraph', installed.packages()[,1])) { stop("Install igraph package.") }

require(igraph)
require(pROC)

source("src/getDCC.R")
source("src/getECC.R")
source("src/getPCC.R")
source("src/getSCC.R")
source("src/getKCC.R")
source("src/rankingByLambda.R")
source("src/getOrderedCummulativeCounts.R")

#' First get the expression matrix that was RMA processed earlier. Used for WDC etc.
#expr.rma.gse3431 <- readRDS(file = "data/Full_RMA_Processed_GSE3431_Expression_Matrix.RDS")
expr.rma <- readRDS(file = "data/Full_RMA_Processed_GSE3431_Expression_Matrix.RDS")

#' Read in PIN graph
hybridG.gse3431 <- readRDS(file = "data/Final_hybrid_network.RDS")

#' Calculate built-in centralities
V(hybridG.gse3431)$degree <- degree(graph = hybridG.gse3431)
V(hybridG.gse3431)$closeness <- closeness(graph = hybridG.gse3431)
V(hybridG.gse3431)$betweenness <- betweenness(graph = hybridG.gse3431)
V(hybridG.gse3431)$eigcent <- evcent(graph = hybridG.gse3431)$vector
V(hybridG.gse3431)$sgc <- subgraph.centrality(graph = hybridG.gse3431, diag = FALSE)

#' Edge betweenness for weighing in centrality
E(hybridG.gse3431)$ebetwn <- edge.betweenness(graph = hybridG.gse3431, directed = F)

#' Edge clustering coefficient (ECC) for weighing in centrality
E(hybridG.gse3431)$ecc <- sapply(X = E(hybridG.gse3431), simplify = T, FUN = function(x){
  getECC(hybridG.gse3431, x)
})

#' Now getPCC for each edge. (Pearson correlation) based on probeset with MAX expression profile sum
E(hybridG.gse3431)$pccMax <- apply(X = igraph::get.edges(hybridG.gse3431, E(hybridG.gse3431)), MARGIN = 1, FUN = function(x){
  getPCC(V(hybridG.gse3431)[x[1]]$Affy.PSs, V(hybridG.gse3431)[x[2]]$Affy.PSs, min = FALSE)
})
#' for missing coexpression due to lack of Affy ids or other issues, will have NAs.
#' We replace the NAs with minimum coexpression value.
E(hybridG.gse3431)$pccMax[ which(is.na(E(hybridG.gse3431)$pccMax)) ] <- min(E(hybridG.gse3431)$pccMax, na.rm = T)

#' Now getPCC for each edge. (Pearson correlation) based on probeset with MIN expression profile sum
E(hybridG.gse3431)$pccMin <- apply(X = igraph::get.edges(hybridG.gse3431, E(hybridG.gse3431)), MARGIN = 1, FUN = function(x){
  getPCC(V(hybridG.gse3431)[x[1]]$Affy.PSs, V(hybridG.gse3431)[x[2]]$Affy.PSs, min = TRUE)
})
#' for missing coexpression due to lack of Affy ids or other issues, will have NAs.
#' We replace the NAs with minimum coexpression value.
E(hybridG.gse3431)$pccMin[ which(is.na(E(hybridG.gse3431)$pccMin)) ] <- min(E(hybridG.gse3431)$pccMin, na.rm = T)

#' Now getDCC for each edges (distance correlation) based on probeset with MAX expression profile sum
E(hybridG.gse3431)$dccMax <- apply(X = igraph::get.edges(hybridG.gse3431, E(hybridG.gse3431)), MARGIN = 1, FUN = function(x){
  getDCC(V(hybridG.gse3431)[x[1]]$Affy.PSs, V(hybridG.gse3431)[x[2]]$Affy.PSs, min = FALSE, index = 1)
})
#' for missing coexpression due to lack of Affy ids or other issues, will have NAs.
#' We replace the NAs with minimum coexpression value.
E(hybridG.gse3431)$dccMax[ which(is.na(E(hybridG.gse3431)$dccMax)) ] <- min(E(hybridG.gse3431)$dccMax, na.rm = T)

#' Now getDCC for each edges (distance correlation) based on probeset with MIN expression profile sum
E(hybridG.gse3431)$dccMin <- apply(X = igraph::get.edges(hybridG.gse3431, E(hybridG.gse3431)), MARGIN = 1, FUN = function(x){
  getDCC(V(hybridG.gse3431)[x[1]]$Affy.PSs, V(hybridG.gse3431)[x[2]]$Affy.PSs, min = TRUE, index = 1)
})
#' for missing coexpression due to lack of Affy ids or other issues, will have NAs.
#' We replace the NAs with minimum coexpression value.
E(hybridG.gse3431)$dccMin[ which(is.na(E(hybridG.gse3431)$dccMin)) ] <- min(E(hybridG.gse3431)$dccMin, na.rm = T)

#' Now getSCC (spearman correlation) for each edge based on probeset with MAX expression profile sum
E(hybridG.gse3431)$sccMax <- apply(X = igraph::get.edges(hybridG.gse3431, E(hybridG.gse3431)), MARGIN = 1, FUN = function(x){
  getSCC(V(hybridG.gse3431)[x[1]]$Affy.PSs, V(hybridG.gse3431)[x[2]]$Affy.PSs, min = FALSE)
})
E(hybridG.gse3431)$sccMax[ which(is.na(E(hybridG.gse3431)$sccMax)) ] <- min(E(hybridG.gse3431)$sccMax, na.rm = T)

#' Now getSCC (spearman correlation) for each edge based on probeset with MIN expression profile sum
E(hybridG.gse3431)$sccMin <- apply(X = igraph::get.edges(hybridG.gse3431, E(hybridG.gse3431)), MARGIN = 1, FUN = function(x){
  getSCC(V(hybridG.gse3431)[x[1]]$Affy.PSs, V(hybridG.gse3431)[x[2]]$Affy.PSs, min = TRUE)
})
E(hybridG.gse3431)$sccMin[ which(is.na(E(hybridG.gse3431)$sccMin)) ] <- min(E(hybridG.gse3431)$sccMin, na.rm = T)

#' Now getKCC (Kendall's tau correlation) for each edge based on probeset with MAX expression profile sum
E(hybridG.gse3431)$kccMax <- apply(X = igraph::get.edges(hybridG.gse3431, E(hybridG.gse3431)), MARGIN = 1, FUN = function(x){
  getKCC(V(hybridG.gse3431)[x[1]]$Affy.PSs, V(hybridG.gse3431)[x[2]]$Affy.PSs, min = FALSE)
})
E(hybridG.gse3431)$kccMax[ which(is.na(E(hybridG.gse3431)$kccMax)) ] <- min(E(hybridG.gse3431)$kccMax, na.rm = T)

#' Now getKCC (spearman correlation) for each edge based on probeset with MIN expression profile sum
E(hybridG.gse3431)$kccMin <- apply(X = igraph::get.edges(hybridG.gse3431, E(hybridG.gse3431)), MARGIN = 1, FUN = function(x){
  getKCC(V(hybridG.gse3431)[x[1]]$Affy.PSs, V(hybridG.gse3431)[x[2]]$Affy.PSs, min = TRUE)
})
E(hybridG.gse3431)$kccMin[ which(is.na(E(hybridG.gse3431)$kccMin)) ] <- min(E(hybridG.gse3431)$kccMin, na.rm = T)


#' Weighted Degree Centrality based on Tang et.al. (2014) "Predicting Essential Proteins Based on Weighted
#' Degree Centrality," IEEE/ACM Trans on Comp Biol and Bioinfo., Vol. 11, No. 2, pp 407-418.
#'
#' $WDC(i) = \sum_{j}^{N_i}{W_{i,j}}$
#' where
#' $W = (ECC \times \lambda) + (PCC \times (1 - \lambda))$
#' and PCC is calculated for probesets with lowest (Min) or highest (Max) expression profile sum.
V(hybridG.gse3431)$wdcMax5 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.5, compOne = "ecc", compTwo = "pccMax")
V(hybridG.gse3431)$wdcMin5 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.5, compOne = "ecc", compTwo = "pccMin")
V(hybridG.gse3431)$wdcMax6 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.6, compOne = "ecc", compTwo = "pccMax")
V(hybridG.gse3431)$wdcMin6 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.6, compOne = "ecc", compTwo = "pccMin")
V(hybridG.gse3431)$wdcMax7 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.7, compOne = "ecc", compTwo = "pccMax")
V(hybridG.gse3431)$wdcMin7 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.7, compOne = "ecc", compTwo = "pccMin")
V(hybridG.gse3431)$wdcMax8 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.8, compOne = "ecc", compTwo = "pccMax")
V(hybridG.gse3431)$wdcMin8 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.8, compOne = "ecc", compTwo = "pccMin")
V(hybridG.gse3431)$wdcMax9 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.9, compOne = "ecc", compTwo = "pccMax")
V(hybridG.gse3431)$wdcMin9 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.9, compOne = "ecc", compTwo = "pccMin")

#' Weighted degree centrality based on Mistry et.al. (2015)
#' $MDC(i) = \sum_{j}^{N_i}{M_{i,j}}$
#' where
#' $M = (ECC \times \lambda) + (DCC \times (1 - \lambda))$
#' and DCC is calculated for probesets with lowest (Min) or highest (Max) expression profile sum.
V(hybridG.gse3431)$mdc1Max5 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.5, compOne = "ecc", compTwo = "dccMax")
V(hybridG.gse3431)$mdc1Min5 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.5, compOne = "ecc", compTwo = "dccMin")
V(hybridG.gse3431)$mdc1Max6 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.6, compOne = "ecc", compTwo = "dccMax")
V(hybridG.gse3431)$mdc1Min6 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.6, compOne = "ecc", compTwo = "dccMin")
V(hybridG.gse3431)$mdc1Max7 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.7, compOne = "ecc", compTwo = "dccMax")
V(hybridG.gse3431)$mdc1Min7 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.7, compOne = "ecc", compTwo = "dccMin")
V(hybridG.gse3431)$mdc1Max8 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.8, compOne = "ecc", compTwo = "dccMax")
V(hybridG.gse3431)$mdc1Min8 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.8, compOne = "ecc", compTwo = "dccMin")
V(hybridG.gse3431)$mdc1Max9 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.9, compOne = "ecc", compTwo = "dccMax")
V(hybridG.gse3431)$mdc1Min9 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.9, compOne = "ecc", compTwo = "dccMin")

#' same as above using ECC and SCC
V(hybridG.gse3431)$mdc2Max5 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.5, compOne = "ecc", compTwo = "sccMax")
V(hybridG.gse3431)$mdc2Min5 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.5, compOne = "ecc", compTwo = "sccMin")
V(hybridG.gse3431)$mdc2Max6 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.6, compOne = "ecc", compTwo = "sccMax")
V(hybridG.gse3431)$mdc2Min6 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.6, compOne = "ecc", compTwo = "sccMin")
V(hybridG.gse3431)$mdc2Max7 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.7, compOne = "ecc", compTwo = "sccMax")
V(hybridG.gse3431)$mdc2Min7 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.7, compOne = "ecc", compTwo = "sccMin")
V(hybridG.gse3431)$mdc2Max8 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.8, compOne = "ecc", compTwo = "sccMax")
V(hybridG.gse3431)$mdc2Min8 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.8, compOne = "ecc", compTwo = "sccMin")
V(hybridG.gse3431)$mdc2Max9 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.9, compOne = "ecc", compTwo = "sccMax")
V(hybridG.gse3431)$mdc2Min9 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.9, compOne = "ecc", compTwo = "sccMin")

#' same as above using ECC and KCC
V(hybridG.gse3431)$mdc3Max5 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.5, compOne = "ecc", compTwo = "kccMax")
V(hybridG.gse3431)$mdc3Min5 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.5, compOne = "ecc", compTwo = "kccMin")
V(hybridG.gse3431)$mdc3Max6 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.6, compOne = "ecc", compTwo = "kccMax")
V(hybridG.gse3431)$mdc3Min6 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.6, compOne = "ecc", compTwo = "kccMin")
V(hybridG.gse3431)$mdc3Max7 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.7, compOne = "ecc", compTwo = "kccMax")
V(hybridG.gse3431)$mdc3Min7 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.7, compOne = "ecc", compTwo = "kccMin")
V(hybridG.gse3431)$mdc3Max8 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.8, compOne = "ecc", compTwo = "kccMax")
V(hybridG.gse3431)$mdc3Min8 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.8, compOne = "ecc", compTwo = "kccMin")
V(hybridG.gse3431)$mdc3Max9 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.9, compOne = "ecc", compTwo = "kccMax")
V(hybridG.gse3431)$mdc3Min9 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.9, compOne = "ecc", compTwo = "kccMin")

#' An alteration to Mistry's weighted centrality
#' $MDC2(i) = \sum_{j}^{N_i}{M_{i,j}}$
#' where
#' $M = (EBC \times \lambda) + (SCC \times (1 - \lambda))$
#' and EBC is edge betweenness centrality, SCC is calculated for probesets
#' with lowest (Min) or highest (Max) expression profile sum.
# V(hybridG.gse3431)$mdc4Max5 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.5, compOne = "ebetwn", compTwo = "sccMax")
# V(hybridG.gse3431)$mdc4Min5 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.5, compOne = "ebetwn", compTwo = "sccMin")
# V(hybridG.gse3431)$mdc4Max6 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.6, compOne = "ebetwn", compTwo = "sccMax")
# V(hybridG.gse3431)$mdc4Min6 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.6, compOne = "ebetwn", compTwo = "sccMin")
# V(hybridG.gse3431)$mdc4Max7 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.7, compOne = "ebetwn", compTwo = "sccMax")
# V(hybridG.gse3431)$mdc4Min7 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.7, compOne = "ebetwn", compTwo = "sccMin")
# V(hybridG.gse3431)$mdc4Max8 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.8, compOne = "ebetwn", compTwo = "sccMax")
# V(hybridG.gse3431)$mdc4Min8 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.8, compOne = "ebetwn", compTwo = "sccMin")
# V(hybridG.gse3431)$mdc4Max8 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.9, compOne = "ebetwn", compTwo = "sccMax")
# V(hybridG.gse3431)$mdc4Min8 <- rankingByLambda(graph = hybridG.gse3431, lambda = 0.9, compOne = "ebetwn", compTwo = "sccMin")

#' Create a data frame to hold all the statistics so that we can do ROC and other tests
graphRankingStats.gse3431 <- data.frame(
                                essential = V(hybridG.gse3431)$Essential,
                                degree = V(hybridG.gse3431)$degree,
                                betweenness = V(hybridG.gse3431)$betweenness,
                                closeness = V(hybridG.gse3431)$closeness,
                                eigcent = V(hybridG.gse3431)$eigcent,
                                sgc = V(hybridG.gse3431)$sgc,

                                wdcMax5 = V(hybridG.gse3431)$wdcMax5,
                                wdcMin5 = V(hybridG.gse3431)$wdcMin5,
                                wdcMax6 = V(hybridG.gse3431)$wdcMax6,
                                wdcMin6 = V(hybridG.gse3431)$wdcMin6,
                                wdcMax7 = V(hybridG.gse3431)$wdcMax7,
                                wdcMin7 = V(hybridG.gse3431)$wdcMin7,
                                wdcMax8 = V(hybridG.gse3431)$wdcMax8,
                                wdcMin8 = V(hybridG.gse3431)$wdcMin8,
                                wdcMax9 = V(hybridG.gse3431)$wdcMax9,
                                wdcMin9 = V(hybridG.gse3431)$wdcMin9,

                                mdc1Max5 = V(hybridG.gse3431)$mdc1Max5,
                                mdc1Min5 = V(hybridG.gse3431)$mdc1Min5,
                                mdc1Max6 = V(hybridG.gse3431)$mdc1Max6,
                                mdc1Min6 = V(hybridG.gse3431)$mdc1Min6,
                                mdc1Max7 = V(hybridG.gse3431)$mdc1Max7,
                                mdc1Min7 = V(hybridG.gse3431)$mdc1Min7,
                                mdc1Max8 = V(hybridG.gse3431)$mdc1Max8,
                                mdc1Min8 = V(hybridG.gse3431)$mdc1Min8,
                                mdc1Max9 = V(hybridG.gse3431)$mdc1Max9,
                                mdc1Min9 = V(hybridG.gse3431)$mdc1Min9,

                                mdc2Max5 = V(hybridG.gse3431)$mdc2Max5,
                                mdc2Min5 = V(hybridG.gse3431)$mdc2Min5,
                                mdc2Max6 = V(hybridG.gse3431)$mdc2Max6,
                                mdc2Min6 = V(hybridG.gse3431)$mdc2Min6,
                                mdc2Max7 = V(hybridG.gse3431)$mdc2Max7,
                                mdc2Min7 = V(hybridG.gse3431)$mdc2Min7,
                                mdc2Max8 = V(hybridG.gse3431)$mdc2Max8,
                                mdc2Min8 = V(hybridG.gse3431)$mdc2Min8,
                                mdc2Max9 = V(hybridG.gse3431)$mdc2Max9,
                                mdc2Min9 = V(hybridG.gse3431)$mdc2Min9,

                                mdc3Max5 = V(hybridG.gse3431)$mdc3Max5,
                                mdc3Min5 = V(hybridG.gse3431)$mdc3Min5,
                                mdc3Max6 = V(hybridG.gse3431)$mdc3Max6,
                                mdc3Min6 = V(hybridG.gse3431)$mdc3Min6,
                                mdc3Max7 = V(hybridG.gse3431)$mdc3Max7,
                                mdc3Min7 = V(hybridG.gse3431)$mdc3Min7,
                                mdc3Max8 = V(hybridG.gse3431)$mdc3Max8,
                                mdc3Min8 = V(hybridG.gse3431)$mdc3Min8,
                                mdc3Max9 = V(hybridG.gse3431)$mdc3Max9,
                                mdc3Min9 = V(hybridG.gse3431)$mdc3Min9,

#                                 mdc4Max5 = V(hybridG.gse3431)$mdc4Max5,
#                                 mdc4Min5 = V(hybridG.gse3431)$mdc4Min5,
#                                 mdc4Max6 = V(hybridG.gse3431)$mdc4Max6,
#                                 mdc4Min6 = V(hybridG.gse3431)$mdc4Min6,
#                                 mdc4Max7 = V(hybridG.gse3431)$mdc4Max7,
#                                 mdc4Min7 = V(hybridG.gse3431)$mdc4Min7,
#                                 mdc4Max8 = V(hybridG.gse3431)$mdc4Max8,
#                                 mdc4Min8 = V(hybridG.gse3431)$mdc4Min8,

                                row.names = V(hybridG.gse3431)$name)

#'
#' Chosen beta, omega based DiffSLc calculations from beta-omega-table.R (results in "doc/AUC_Table.xlsx")
#'
om <- 0.9
graphRankingStats.gse3431$NT1_diffslc <- (om*V(hybridG.gse3431)$eigcent) + ((1-om)*graphRankingStats.gse3431$wdcMin8)

om <- 0.3
graphRankingStats.gse3431$NT2_diffslc <- (om*V(hybridG.gse3431)$eigcent) + ((1-om)*graphRankingStats.gse3431$mdc1Min8)

om <- 0.9
graphRankingStats.gse3431$NT3_diffslc <- (om*V(hybridG.gse3431)$eigcent) + ((1-om)*graphRankingStats.gse3431$mdc2Min7)

#'
#' Chosen beta, omega based on DiffSLc calc from beta-omega-table.R (results in "doc/AUC_Table.xlsx")
#'
om <- 0.9
graphRankingStats.gse3431$NF1_diffslc <- (om*V(hybridG.gse3431)$eigcent) + ((1-om)*graphRankingStats.gse3431$wdcMin8)

om <- 0.3
graphRankingStats.gse3431$NF2_diffslc <- (om*V(hybridG.gse3431)$eigcent) + ((1-om)*graphRankingStats.gse3431$mdc1Min8)

om <- 0.9
graphRankingStats.gse3431$NF3_diffslc <- (om*V(hybridG.gse3431)$eigcent) + ((1-om)*graphRankingStats.gse3431$mdc2Min8)

# Get AUC and ROC.
ROCs <- sapply(X = colnames(graphRankingStats.gse3431), simplify = T, FUN = function(aCol, df, countColumn){
  if(aCol != countColumn) {
    roc(df[,countColumn] ~ df[,aCol], auc = TRUE, ci = TRUE)
  } else {
    NULL
  }
}, countColumn = "essential", df = graphRankingStats.gse3431)
