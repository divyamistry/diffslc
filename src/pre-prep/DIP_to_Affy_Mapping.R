#' @author Divya Mistry
#' @description Fill in probeset ids for DIP interactors
#' 

#' Read Yeast annotation
YG_S98Annot <- read.csv(file = "data/YG_S98.na34.annot.csv", header = T, skip=19, na.string = "---", stringsAsFactors = F)
#' Read DIP-to-Uniprot fetched from Uniprot/SwissProt KB
DIPtoUniprotMapping <- read.delim(file = "data/DIP-to-Uniprot-Matched-IDs.tab", header = T, na.string="", stringsAsFactors = F)
colnames(DIPtoUniprotMapping) <- c("UniProtKB", "DIP", "Entry.name", "Protein.names", "Gene.names", 
                                   "Length", "Gene.names.ORF", "Gene.names.ordered.locus", 
                                   "Gene.names.synonym", "Gene.names.primary", "Cross.ref.RefSeq", 
                                   "Cross.ref.EnsembleFungi", "Cross.ref.GeneID", "Cross.ref.CYGD", "Cross.ref.SGD")

#' DipMapper() maps a DIP id to up to 4 affy probesets.
source('src/pre-prep/DipMapper.R')
DIPtoUniprotMapping$Probe.Set.IDs <- sapply(X = DIPtoUniprotMapping$DIP, simplify = T, USE.NAMES = T, FUN = function(x) DipMapper(x))

#' Write the processed data frame for future use.
write.table(x = DIPtoUniprotMapping, file = "data/DIPtoAffy_with_additionalAnnotations.tsv", quote = F, sep = "\t", na = "---", row.names = F, col.names = T)
