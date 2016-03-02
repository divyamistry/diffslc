#' @author Divya Mistry
#' @description Generate the MDC ranking using given lambda value
#' 
#' @param graph igraph for which ranking is to be calculated
#' @param lambda the weighing factor in [0,1] range
#' @param compOne string, the first component of the weighing formula (See Details). Pre-calculated as edge attribute.
#' @param compTwo string, the second component of the weighing formula (See Details). Pre-calculated as edge attribute.
#' 
#' @return centrality ranking score in the same order as V(graph)
#' 
#' @details MDC is calculated using the following
#' Weighted degree centrality based on Mistry et.al. (2015)
#' $MDC(i) = \sum_{j}^{N_i}{M_{i,j}}$
#' where
#' $M = (compOne \times \lambda) + (compTwo \times (1 - \lambda))$
#' compOne and compTwo are edge weighing values from a node's incident edges. For example, compOne
#' and compTwo could be DCC and ECC of edges.
#' 
rankingByLambda <- function(graph, lambda, compOne, compTwo) {
  # check if given components/attributes actually exist
  if(!(compOne %in% list.edge.attributes(graph))) {
    stop(paste(compOne,"isn't one of the edge attributes of graph"))
  }
  if(!(compTwo %in% list.edge.attributes(graph))) {
    stop(paste(compTwo, "isn't one of the edge attributes of graph"))
  }
  # return the score
  sapply(X = V(graph), simplify = TRUE, FUN = function(v){
    incidentEdges <- E(graph)[unlist(graph[[v,,edges=TRUE]])]
    sum(
      (lambda * get.edge.attribute(graph = graph, name = compOne, index = incidentEdges))
      +
      ((1-lambda) * get.edge.attribute(graph = graph, name = compTwo, index = incidentEdges))
    )
  })
}
# rm(graph,lambda,compOne,compTwo)

# incidentEdges <- E(hybridG)[unlist(hybridG[[2,,edges=TRUE]])]

rankingByAverage <- function(graph, compOne) {
  # check if given components/attributes actually exist
  if(!(compOne %in% list.edge.attributes(graph))) {
    stop(paste(compOne,"isn't one of the edge attributes of graph"))
  }
  # return the score
  sapply(X = V(graph), simplify = TRUE, FUN = function(v){
    incidentEdges <- E(graph)[unlist(graph[[v,,edges=TRUE]])]
    length(incidentEdges) * sum(get.edge.attribute(graph = graph, name = compOne, index = incidentEdges))
  })
}