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

#' Data frame to plot NT1,2,3
measuresToCalc <- c("NT1_diffslc", "NT2_diffslc", "NT3_diffslc") #pcor, scor, dcor

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
  labs(title = "ROC curves for NT1, NT2, NT3 networks", x = "FPR", y = "TPR") +
  geom_line(size = 0.5) +
  scale_color_discrete(name = "Centrality",
                       breaks = levels(resultDF$Centrality),
                       labels = paste(c("NT1", "NT2", "NT3"),", AUC: ",format(resultAUC,digits = 4),sep =''))
plot(myplot)

#' plot %improvement from DC in Top-1, 5, 10, 15, 20, 25% of matches.
# Calculate cummulative counts of matched TRUE values in "essential" column
# and get difference from degree counts to get %change values
cummulCounts <- getOrderedCummulativeCounts(df = graphRankingStats, countColumn = "essential", decreasing = T)

# plot cummulative counts
pctMatchForPlotting <- c(0.01, 0.05, 0.1, 0.15, 0.20, 0.25) * vcount(hybridG) # length(which(graphRankingStats$essential))
for(pct in pctMatchForPlotting) {
  degFixPt <- cummulCounts[pct,"degree"] # count of matched essential genes using degree centrality
  yminVal <- min(cummulCounts[pct,c("degree",measuresToCalc)])
  ymaxVal <- max(cummulCounts[pct,c("degree",measuresToCalc)])
  improvementPcts <- cummulCounts[pct,measuresToCalc] - cummulCounts[pct,"degree"]
  # improvementPcts <- abs(cummulCounts[pct,measuresToCalc] - cummulCounts[pct,"degree"])/cummulCounts[pct,"degree"] # % improvement from Degree centrality
  barplot(height = as.matrix(improvementPcts),
          legend.text = F, horiz = F, ylim = c(0,50),
          main = paste("Matched counts at ", 100*pct/vcount(hybridG), "%", sep=''), las = 2,
          ylab = "Additional essential genes detected by each method", names.arg = c("NT1", "NT2", "NT3"))
}

#' All tens of measures compared together as a stacked barchart
# barplot(height = as.matrix(cummulCounts[pctMatchForPlotting,c("degree",measuresToCalc)]),
#         legend.text = F, horiz = F, ylim = c(0,2000),
#         main = "Matched counts at 1-25%", offset = 0, las=2, 
#         names.arg = colnames(graphRankingStats)[-1]) #-1 b/c no need to plot "essential" column

