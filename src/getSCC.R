#' @author Divya Mistry
#' @description Function to get spearman correlation based on DIP interactors
#' 
#' @param nodeApsets probesets for one of the incident nodes for edge of interest
#' @param nodeBpsets probesets for the other incident node to edge of interest
#' @param min logical, TRUE for using correlation between probesets with lowest expression profile sum, and
#'                     FALSE (default) for using correlation between probesets with highest expression profile sum.
#' 
#' @details Requires the expression data loaded in current environment
#' @return pairwise pearson correlation value for most abundant probesets corresponding to the DIP interactors represted by nodeA--nodeB edge.
#' 
message("NOTE: getSCC() requires expression matrix in the current environment")

# nodeApsets<-V(connGraph)[287]$Affy.PSs; nodeBpsets<-V(connGraph)[288]$Affy.PSs
getSCC <- function(nodeApsets, nodeBpsets, min = FALSE) {
  # if either of the probeset lists are NA, return NA. The NAs will processed for least coexpr
  #   value at the end of all coexprEdgeWt
  if(is.na(nodeApsets) || is.na(nodeBpsets)) {
    NA
  } else {
    # split the nodeApsets
    splitPs <- unlist(strsplit(x = nodeApsets, split = ";", fixed = T))
    # get the probeset with least expression amount of the group
    if(length(splitPs) == 1) {
      splitA <- splitPs
    } else {
      leastExpr <- rowSums(x = expr.rma[splitPs,])
      if(min == TRUE) {
        splitA <- splitPs[which(leastExpr == min(leastExpr))[1]] # index by [1] to pick first one if multiple expressions add up to minimum
      } else {
        splitA <- splitPs[which(leastExpr == max(leastExpr))[1]] # index by [1] to pick first one if multiple expressions add up to minimum
      }
    }
    
    # split the nodeBpsets
    splitPs <- unlist(strsplit(x = nodeBpsets, split = ";", fixed = T))
    # get the probeset with least expression amount of the group
    if(length(splitPs) == 1) {
      splitB <- splitPs
    } else {
      leastExpr <- rowSums(x = expr.rma[splitPs,])
      if(min == TRUE) {
        splitB <- splitPs[which(leastExpr == min(leastExpr))[1]] # index by [1] to pick first one if multiple expressions add up to minimum
      } else {
        splitB <- splitPs[which(leastExpr == max(leastExpr))[1]] # index by [1] to pick first one if multiple expressions add up to minimum
      }
    }
    
    # return pearson correlation for the highest expressed probesets
    cor(x = expr.rma[splitA,], y = expr.rma[splitB,], method = "spearman")
  }
}