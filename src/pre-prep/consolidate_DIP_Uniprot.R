#' @author Divya Mistry
#' @description consolidate and trim the interaction dataset to match the data from UniProt
#' 

#' I got DIP->Uniprot mapped yeast reads from uniprot website. I read and hold that data
uniprotMappedDip <- read.delim("data/DIP-to-Uniprot-Matched-IDs.tab", na.string = '', stringsAsFactors = F)

#' DIP interactions that have uniprot known
# mappedDIP <- intrxns$ID.interactor.A[intrxns$ID.interactor.A %in% uniprotMappedDip$yourlist.M2015040895A3I6Y1GK]
mappedDIPViaUniprotA <- intrxns$UniprotKB.interactor.A %in% uniprotMappedDip$Entry
mappedDIPViaUniprotB <- intrxns$UniprotKB.interactor.B %in% uniprotMappedDip$Entry
length(which(mappedDIPViaUniprotA == mappedDIPViaUniprotB)) # 22012 elements match

#' columns of YGS98 annotation data that have filled in info
naCountInYG <- apply(X = YG_S98Annot, MARGIN = 2, FUN = function(aCol){ length(which(is.na(aCol)))})
YG_S98Annot_mod <- YG_S98Annot[,which(naCountInYG == 0)]
YG_S98Annot_mod <- YG_S98Annot_mod[,-c(2:6,8,10:12)] # useless info for now
#' remove control probes from YG_S98, There are 60 of them, all with AFFX- prefix in Representative ID
# length(unique(grep(pattern = "AFFX-", x = YG_S98Annot$Probe.Set.ID, fixed = T))) # gives 60, which is what we need
YG_S98Annot_mod <- YG_S98Annot_mod[-grep(pattern = "AFFX-", x = YG_S98Annot$Probe.Set.ID, fixed = T),]

#' Read in essential genes list from DEG
fullDEG <- read.delim(file = "data/degannotation-e.dat", stringsAsFactors = F, na.string = "-")
#' filter yeast form Essential genes list
fullDEG <- subset(x = fullDEG, X.Organism == "Saccharomyces cerevisiae")

delme3<-sapply(X = YG_S98Annot_mod$Representative.Public.ID, simplify = T, FUN = function(rpi){ grep(pattern = rpi, x = uniprotMappedDip$Gene.names, fixed = T) })

#' which DIP can be mapped to probesets
#' 
delme <- sapply(X = uniprotMappedDip$DIP, simplify = T, FUN = function(DIPid){
  idx <- which(uniprotMappedDip$DIP == DIPid) #multiple results?
  upID <- uniprotMappedDip$Entry[idx]
  ys_upID_idxs <- grep(pattern = upID, x = YG_S98Annot$SwissProt, fixed = TRUE) #multiple results?
  psIDs <- as.character(YG_S98Annot$Probe.Set.ID[ys_upID_idxs])
})
delme2 <- sapply(X = delme, simplify = T, FUN = function(dm){ length(dm) })
length(which(delme2 != 0)) # 4609 elements.
#' We started with 5077 DIP based nodes in our graph, and are only able to find
#' 4609 probeset ids to match with that. The other ~9% of the DIP id nodes will
#' have to suffer from min value choosing
