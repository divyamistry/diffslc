#' b is for beta parameter
#' o is for omega parameter
#' 

omegas <- seq(from = 0.1, to = 0.9, by = 0.1)
# beta's [0.5,0.8] have already been computed earlier. I'll just pull that from precomputed results

#' for pCor (b,o) -> (0.8,0.9)
cols <- colnames(graphRankingStats)[7:16]
pcorAUCmatrix<-matrix(nrow = length(omegas), ncol = length(cols), dimnames = list(as.character(omegas), cols))
for(b in cols){
  for(o in omegas){
    diffslc <- (o*V(hybridG)$eigcent) + ((1-o)*graphRankingStats[,b])
    pcorAUCmatrix[as.character(o),as.character(b)] <- roc(graphRankingStats$essential ~ diffslc)$auc
  }
}

#' for sCor (b,o) -> (0.8,0.5)
cols <- colnames(graphRankingStats)[17:26]
scorAUCmatrix<-matrix(nrow = length(omegas), ncol = length(cols), dimnames = list(as.character(omegas), cols))
for(b in cols) {
  for(o in omegas){
    diffslc <- (o*V(hybridG)$eigcent) + ((1-o)*graphRankingStats[,b])
    scorAUCmatrix[as.character(o),as.character(b)] <- roc(graphRankingStats$essential ~ diffslc)$auc
  }
}

#' for kCor (b,o) -> (0.7,0.9) or (0.8,0.9)
cols <- colnames(graphRankingStats)[27:36]
kcorAUCmatrix<-matrix(nrow = length(omegas), ncol = length(cols), dimnames = list(as.character(omegas), cols))
for(b in cols){
  for(o in omegas){
    diffslc <- (o*V(hybridG)$eigcent) + ((1-o)*graphRankingStats[,b])
    kcorAUCmatrix[as.character(o),as.character(b)] <- roc(graphRankingStats$essential ~ diffslc)$auc
  }
}

#' for dCor (b,o) -> (0.8,0.9) or (0.8,0.8)
cols <- colnames(graphRankingStats)[37:46]
dcorAUCmatrix<-matrix(nrow = length(omegas), ncol = length(cols), dimnames = list(as.character(omegas), cols))
for(b in cols){
  for(o in omegas){
    diffslc <- (o*V(hybridG)$eigcent) + ((1-o)*graphRankingStats[,b])
    dcorAUCmatrix[as.character(o),as.character(b)] <- roc(graphRankingStats$essential ~ diffslc)$auc
  }
}

