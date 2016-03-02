#' @author Divya Mistry
#' @description Function to scale the given vector in [0,1] range.
#'              For a given vector X, This is achieved by $(x_i - min(X)) / (max(X) - min(X))$
#'              If $min(X) == max(X)$ is TRUE, function returns 1.
#' @param X numeric vector to rescale
#' @return X rescaled in [0,1]. 
scaleFromZeroToOne <- function(X) {
  if(min(X) == max(X)) {
    rep(1,length(X))
  } else {
    (X - min(X)) / (max(X) - min(X))
  }
}