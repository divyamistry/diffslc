#' @author Divya Mistry
#' @description Read the DIP (Database of Interacting Proteins http://dip.doe-mbi.ucla.edu/dip/Main.cgi)
#'              yeast data from 2015/01/02 update. Create a resulting minimal 
#'              dataframe with interaction data
#'              

# dip_data will have 22838 interactions to start
dip_data <- read.table(file = "data/yeastDIP_Scere20150101.txt", header = TRUE,
                       sep = "\t", na.string = "-", stringsAsFactors = FALSE, strip.white = TRUE)

#' for whatever reason, the PSI-Mi Tab 2.5 formatted file from DIP has last two columns 
#' with empty or NA data. And the columns don't have names, so it's not a well-formatted
#' tsv file. I've manually added two column names in the raw text file, and now I'm going
#' to just discard that data
#' 
dip_data <- dip_data[, -c(3:6,8,17,18)]

#' Now remove interactions involving non-yeast species. dip_data will ahve 22680 interactions after this.
yeast_intrxns <- intersect(grep("Saccharomyces cerevisiae", dip_data$Taxid.interactor.A), 
                           grep("Saccharomyces cerevisiae", dip_data$Taxid.interactor.B))
dip_data <- dip_data[yeast_intrxns,]

#' Grab the refseq and uniprotkb ids so that I can compare them against YG_S98 annotations and get
#' corresponding probesets for coexpression values. Just a personal opinion that MI-TAB is not
#' the nicest format to work with. For a delimiter separated value file, this is quite painful
#' to have to keep splitting values. 
#' 
Interactor.A<-strsplit(x = dip_data$ID.interactor.A, split = "[|:]", fixed = FALSE)
Interactor.B<-strsplit(x = dip_data$ID.interactor.B, split = "[|:]", fixed = FALSE)

splittedA <- sapply(X = Interactor.A, simplify = T, FUN = function(yy){
  if(length(yy)==1) { # only DIP exists
    c(yy[1],NA,NA) #return NA for missing values
  } else if (length(yy)==3) { # either refseq or uniprotkb is missing
    if(yy[2] == "refseq") {
      c(yy[1],yy[3],NA)
    } else {
      c(yy[1],NA,yy[3])
    }
  } else if (length(yy)==5) { # all three DIP, refseq, uniprotkb available
    c(yy[1],yy[3],yy[5])
  } else { # there's more info than expected. report the result
    warning(paste("Couldn't correctly parse this row: ",yy))
  }
})
splittedA <- t(splittedA)

splittedB <- sapply(X = Interactor.B, simplify = T, FUN = function(yy){
  if(length(yy)==1) { # only DIP exists
    c(yy[1],NA,NA) #return NA for missing values
  } else if (length(yy)==3) { # either refseq or uniprotkb is missing
    if(yy[2] == "refseq") {
      c(yy[1],yy[3],NA)
    } else {
      c(yy[1],NA,yy[3])
    }
  } else if (length(yy)==5) { # all three DIP, refseq, uniprotkb available
    c(yy[1],yy[3],yy[5])
  } else { # there's more info than expected. report the result
    warning(paste("Couldn't correctly parse this row: ",yy))
  }
})
splittedB <- t(splittedB)

#' Interactions for which refseq or uniprot references don't exist.
#' We need either to match the data against YG_S98 annotations 
#' 
View(dip_data[ which(!is.na(splittedA[,2]) | !is.na(splittedB[,2])), ]) #refseq
View(dip_data[ which(!is.na(splittedA[,3]) | !is.na(splittedB[,3])), ]) #uniprot

#' Create a data frame and save the interaction info in file
#' 
intrxns <- cbind(splittedA,splittedB)
colnames(intrxns) <- c("ID.interactor.A","RefSeq.interactor.A","UniprotKB.interactor.A",
                       "ID.interactor.B","RefSeq.interactor.B","UniprotKB.interactor.B")
write.table(x = intrxns, file = "data/minimal_DIP_interaction_data.csv", quote = F, sep = ",", na = "-", row.names = FALSE, col.names = TRUE)
#' command to read the file later
#' intrxns2 <- read.table(file = "data/minimal_DIP_interaction_data.csv", header = T, sep = ",", na.strings = "-", stringsAsFactors = F)

#' remove files created in this file
rm(dip_data, intrxns, splittedA, splittedB, Interactor.A, Interactor.B, yeast_intrxns)
