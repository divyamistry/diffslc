#' @name getOrderedCummulativeCounts
#' @description Function to generate cummulative counts of TRUE values in a given column of
#'              data frame, based on other sorted columns.
#' 
#' @param df data frame to count sorted results in.
#' @param countColumn column in which the number of TRUE values will be counted.
#' @param decreasing logical, TRUE (default) for counts based on descending order of columns in 
#'                   given data frame. FALSE for ascending ordered values to be counted against
#'                   given countColumn.
#' @return data frame of cummulative running counts of each of the columns for the number of TRUE
#'         values in given countColumn
getOrderedCummulativeCounts <- function(df = graphRankingStats, countColumn = "essential", decreasing = "T") {
  calculatedDf <- data.frame(row.names = 1:length(which(df[,countColumn]))) # prepare a data frame to add counts to
  for(i in 1:ncol(df)) { # iterate through the columns
    if(colnames(df)[i] != countColumn){ # skip the column that we're counting against.
      # order by current column, and get corresponding countColumn layout. We'll be cummulatively counting that
      dfOrdered <- df[ order(df[,i], decreasing = decreasing), countColumn ]
      # go through the rows to start adding up counts of TRUE from the countColumn, and attach them to the result dataframe
      calculatedDf <- cbind(calculatedDf, sapply(X = 1:length(which(dfOrdered)), simplify = T,
           FUN = function(countPos, dfOrdered){ # countPos is the current position in the counts of countColumn
             length( which( dfOrdered[1:countPos]))
           },
           dfOrdered
        )
      )
      # let's fix the column name to the original column name of given data frame
      colnames(calculatedDf)[ncol(calculatedDf)]<-colnames(df)[i]
    }
  }
  # return the resulting data frame.
  calculatedDf
}