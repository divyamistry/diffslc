#' @author Divya Mistry
#' @description This file is used to read in GSE3431 expression CEL files, normalize the data
#' and then save the processed data so that it can be read in later to get
#' pairwise coexpression data
#' 

require(affy)
require(microbenchmark)

#' I downloaded the raw CEL files from author's site http://moment.utmb.edu/cgi-bin/dload.cgi

#' I read in CEL files so that I can RMA normalize
gse3431 <- ReadAffy(celfile.path = "data/GSE3431/")
rma.gse <- rma(gse3431)

#' expression matrix
expr.rma <- exprs(rma.gse)

#' sanity-check for normalized expressions
boxplot(expr.rma, las=2)

saveRDS(object = expr.rma, file = "data/Full_RMA_Processed_GSE3431_Expression_Matrix.RDS", compress = TRUE)

#' According to Li et.al. (2014) and Tang et.al. (2014) the data has 6777 "genes"
#' This is a result of noticing that multiple probesets seem to map to the same
#' transcribed region. A representative ID for this is obtained from reading
#' Affymetrix YG_S98 annotation file. The file has been downloaded in "data/GSE3431/"
#' In the following lines, I read the annotation file, and verify that there are
#' 6777 transcribed region, i.e. "Representative Public ID" as described at
#' http://www.affymetrix.com/support/technical/manual/taf_manual.affx#probe_design_information
#' 
YG_S98Annot <- read.csv(file = "data/YG_S98.na34.annot.csv", header = T, skip=19, na.string = "---")
length(levels(YG_S98Annot$Representative.Public.ID))

#' Let's make sure probeset names from annotation file and experiment data are same
#' Following test makes sure that not only are the probesets the same, but that 
#' they are in the same order. This comes in handy when referencing things by index
#' between different data sources.
#' 
which(YG_S98Annot$Probe.Set.ID != rownames(rma.gse))

#' Now we read the expression profiles of probesets that map to same transcribed region.
#' We keep the probeset representing the highest mean expression over all 36 conditions.
#' This is simply to avoid breaking apart expressions that probably belongs to
#' one gene, into multiple "genes". This is not a perfect solution, but it's good enough.
#' 
transcribedRegions <- levels(YG_S98Annot$Representative.Public.ID)
keeperPS <- sapply(X = transcribedRegions, simplify = T, FUN = function(x){ # took around a minute and half to run
  currResult <- YG_S98Annot[ YG_S98Annot$Representative.Public.ID == x, ]
  if(dim(currResult)[1]>1) { # these are the probesets with multiple choices
    # the probesets with highest total expression will end up having highest average
    rSums <- rowSums(expr.rma[ currResult$Probe.Set.ID, ])
    # we save the highest contributor probeset
    psKeeper <- which(rSums == max(rSums))
    currResult$Probe.Set.ID[psKeeper]
  } else {
    # these are the unique choices, so we return them as-is
    currResult$Probe.Set.ID
  }
})

#' Save the list of unique transcribed regions separately for more use.
saveRDS(object = keeperPS, file = "data/gse3431_uniqueProbeSets_represented_by_NetAffx.RDS", ascii = T)

#' remove unnecessary objects
rm(expr.rma, rma.gse, gse3431, keeperPS, transcribedRegions, YG_S98Annot)
