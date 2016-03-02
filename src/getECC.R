#' @author Divya Mistry
#' @description Function to calculate Edge clustering coefficient for a given edge as
#'              defined in Wang et.al. "Identification of Essential Proteins Based on 
#'              Edge Clustering Coefficient," IEEE/ACM Trans. Computational Biology 
#'              and Bioinformatics, vol. 9, no. 4, pp. 1070-1080, July/Aug. 2012
#' 


#' @name getECC
#' @description calculate ECC using the following formula
#'              $ECC(u,v) = \frac{| N_u \subset N_v | + 1}{min\{d_u, d_v\}}$
#'              where u and v are vertices incident to edge. N_u and N_v are 
#'              neighbors of nodes u and v. d_u and d_v are degrees of the nodes.
#' @param graph - an igraph object for given nodes
#' @param e - edge for which the ECC is to be computed
#' @return a scalar, Edge clustering coefficient
getECC <- function(graph, e) {
  ve <- get.edge(graph, e)
  numerator <- length(intersect(neighbors(graph, v = ve[1], mode = 1), neighbors(graph, v = ve[2], mode = 1))) + 1
  denominator <- min(degree(graph = graph, v = ve[1], loops = F), degree(graph = graph, v = ve[2], loops = F))
  ECC <- numerator / denominator
  ECC
}