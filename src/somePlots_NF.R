# ROC plots aes(specificities, sensitivities)
library(ggplot2)
library(ROCR) #for precision recall curves.

centNames<-colnames(graphRankingStats);
for(i in 1:length(ROCs)) ROCs[[i]]$centName <- centNames[i]

ROCdf<-lapply(X = ROCs, FUN = function(r){
  if(!is.null(r$call)) {
    data.frame(FPR = (1 - r$specificities), TPR = r$sensitivities, Centrality = r$centName)
  } else {
    NULL
  }
})

#' ROC calls for post-processed data.

#' Data and process to plot NF1,2,3
#' NF1,2,3 ROC curve figure
measuresToCalc <- c("NF1_diffslc", "NF2_diffslc", "NF3_diffslc") #pcor, scor, dcor

resultDF <- data.frame()
for(i in measuresToCalc) {
  if(!is.null(ROCdf[[i]])) {
    resultDF <- rbind(resultDF, ROCdf[[i]])
  }
}

resultAUC <- c()
for(i in measuresToCalc){
  if(!is.null(ROCs[[i]]$call)) {
    resultAUC <- c(resultAUC, ROCs[[i]]$auc)
  }
}
names(resultAUC) <- levels(resultDF$Centrality)

myplot <- ggplot(data = resultDF, aes(FPR, TPR, color=Centrality)) + 
  labs(title = "ROC curves for NF1, NF2, NF3 networks", x = "FPR", y = "TPR") +
  geom_line(size = 0.5) +
  scale_color_discrete(name = "Centrality",
                       breaks = levels(resultDF$Centrality),
                       labels = paste(c("NF1", "NF2", "NF3"),", AUC: ",format(resultAUC,digits = 4),sep =''))
plot(myplot)

#' Plot improvements in for NF1,2,3 for top-1, 5, 10, 15, 20, 25% 
# pctMatchForPlotting <- c(0.01, 0.05, 0.1, 0.15, 0.20, 0.25) * vcount(hybridG) # length(which(graphRankingStats$essential))
for(pct in pctMatchForPlotting) {
  degFixPt <- cummulCounts[pct,"degree"] # count of matched essential genes using degree centrality
  yminVal <- min(cummulCounts[pct,c("degree",measuresToCalc)])
  ymaxVal <- max(cummulCounts[pct,c("degree",measuresToCalc)])
  improvementPcts <- cummulCounts[pct,measuresToCalc] - cummulCounts[pct,"degree"]
  # improvementPcts <- abs(cummulCounts[pct,measuresToCalc] - cummulCounts[pct,"degree"])/cummulCounts[pct,"degree"] # % improvement from Degree centrality
  barplot(height = as.matrix(improvementPcts),
          legend.text = F, horiz = F, ylim = c(0,50),
          main = paste("Matched counts at ", 100*pct/vcount(hybridG), "%", sep=''), las = 2,
          ylab = "Additional essential genes detected by each method", names.arg = c("NF1", "NF2", "NF3"))
}
