#' @author Divya Mistry
#' @description Find out the discrepancies between DIP and YG_S98, and create a dataframe
#'              with all the available mappings. Save the mappings for future use
#' 

#' Read in the minimal interaction information prepared from full DIP data
#' 
intrxns <- read.table(file = "data/minimal_DIP_interaction_data.csv", header = T, sep = ",", na.strings = "-", stringsAsFactors = F)
splittedA <- intrxns[,1:3]
splittedB <- intrxns[,4:6]

#' some of the SwissProt column values in YG_S98 annotation have multiple listings
#' which requries a substring comparison instead of full string comparison.
#' We can split by the delimiter '///' and compile full list of available values and 
#' then see which refseq or uniprot/swissprot id's are missing
#' Let's see what we can't to find in YG_S98
#' 
splitSwissProt <- strsplit(x = as.character(YG_S98Annot$SwissProt), split = " /// ", fixed = T)
splitSwissProt <- unlist(splitSwissProt)
length(which(as.character(splittedA[,3]) %nin% as.character(splitSwissProt))) # 515 elements
length(which(as.character(splittedB[,3]) %nin% as.character(splitSwissProt))) # 555 elements

splitRSprot <- strsplit(x = as.character(YG_S98Annot$RefSeq.Protein.ID), split = " /// ", fixed = T)
splitRSprot <- unlist(splitRSprot)
length(which(as.character(splittedA[,2]) %nin% as.character(splitRSprot))) # 183 elements
length(which(as.character(splittedB[,2]) %nin% as.character(splitRSprot))) # 226 elemtns

splitRStran <- strsplit(x = as.character(YG_S98Annot$RefSeq.Transcript.ID), split = " /// ", fixed = T)
splitRStran <- unlist(splitRStran)
length(which(as.character(splittedA[,2]) %nin% as.character(splitRStran))) # 22465 elements
length(which(as.character(splittedB[,2]) %nin% as.character(splitRStran))) # 22525 elemtns

#' Let's start removing self-loops since we don't need them for this analysis
#' Also, there don't seem to be repeated interactions based on interactor ID's.
intrxns <- intrxns[-which(intrxns$ID.interactor.A == intrxns$ID.interactor.B) ,]




