#' @author Divya Mistry
#' @description Function to map DIP interactor ID to corresponding up to 4 YG_S98 probeset ids

#' @name DipMapper
#' @description Function to map DIP interactor ID to corresponding up to 4 YG_S98 probeset ids
#' 
#' @param DIPstr DIP interactor ID for yeast
#' 
#' @details Expecting following two data farmes already available in current environment
#'          - YG_S98Annot YG_S98 Release 34 annotation parsed in a data frame with stringsAsFactors=FALSE.
#'          - DIPtoUniprotMapping data frame with DIP-UniProt-OrderedLocus mappings. Available in data/
#' 
#' @return a semi-colon separated list of up to four YG_S98 probeset ids mapped to given DIP interactor. 
#'         For unmatched or more than 4 matches, NA is returned.
#'         
message("NOTE: DipMapper expects YG_S98Annot and DIPtoUniprotMapping data frames availabe in current environment.")
DipMapper <- function(DIPstr) {
  orderedLocus <- DIPtoUniprotMapping$Gene.names.ordered.locus[ which(DIPtoUniprotMapping$DIP == DIPstr) ]
  if( (length(orderedLocus)==0) || is.na(orderedLocus) ) { #if no ordered locus. Throw warning and return NA
    warning(paste("DIP didn't match against any known ordered.locus for", DIPstr))
    NA
  } else {
    splitOL <- unlist(strsplit(x = orderedLocus, split = "; ", fixed = T))
    matchedIdx <- unlist(sapply(X = splitOL, simplify = T, USE.NAMES = T, FUN = function(ol){
      matchedLoc <- grep(ol, YG_S98Annot$Representative.Public.ID, ignore.case = T)
      if(length(matchedLoc) == 0) { # if Representive Public ID wasn't found, look in Ensembl
        matchedLoc <- grep(ol, YG_S98Annot$Ensembl, ignore.case = T)
      }
      matchedLoc
    }))
    # if less than 5 matchedIdx's were found, we return all probesets 
    # otherwise we return NA. The limit of 5 was determined based on matches
    # against the Representative.Public.ID. All the ordered.locus that match
    # exactly to Representative Public ID, end up having at most 4 matches.
    # Therefore, I limit complexity to the same number of matches for Ensembl as well.
    matchedIdx <- unique(matchedIdx[!is.na(matchedIdx)])
    if((length(matchedIdx) == 0) || (length(matchedIdx) > 4)) {
      NA
    } else {
      matchedProbesets <- YG_S98Annot$Probe.Set.ID[ matchedIdx ]
      paste(matchedProbesets, collapse=';')
    }
  }
}

#' @name DipEssential
#' @description Map DIP interactor name string to essentiality as 
#' @param DIPstr DIP interactor ID
#' @return TRUE if DIP interactor is essential, FALSE otherwise
#' 
DipEssential <- function(DIPstr){
  queryGenes <- DIPtoUniprotMapping$Gene.names[ which(DIPtoUniprotMapping$DIP == DIPstr) ]
  if( (length(queryGenes) == 0) || is.na(queryGenes) ) { #if no Gene.names were found for DIP
    warning(paste("DIP didn't match against any known Gene.names for ", DIPstr, ". DIP will be considered non-essential by default.", sep = ""))
    FALSE
  } else {
    splitGenes <- unlist(strsplit(queryGenes, split= " "))
    matchedEssential <- sapply(X = splitGenes, simplify = T, FUN = function(x){ # for each of the query list
      unlist( # unlist to easily count number of essential gene matches
        sapply(X = x, simplify = T, FUN = function(y){ # for each gene name within list item
          grep(y, fullDEG$X.Gene_Name) # if query gene name matches DEG gene names
        })
      )
    })
    #which(sapply(X = matchedEssential, simplify = T, FUN = length) > 0)
    if( length(unlist(matchedEssential)) > 0 ){ # essential gene match was found
      TRUE
    } else { # essential gene match wasn't found
      FALSE
    }
  }
}

