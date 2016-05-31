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

#' data frame to plot N0 centrality stats
measuresToCalc <- c("degree","betweenness","closeness","eigcent","sgc")
resultDF <- data.frame()
for(i in measuresToCalc) {
  if(!is.null(ROCdf[[i]])) {
    resultDF <- rbind(resultDF, ROCdf[[i]])
  }
}

resultAUC <- c()
for(i in measuresToCalc) {
  if(!is.null(ROCs[[i]]$call)) {
    resultAUC <- c(resultAUC, ROCs[[i]]$auc)
  }
}
names(resultAUC)<-levels(resultDF$Centrality)

myplot <- ggplot(data = resultDF, aes(FPR, TPR, color=Centrality)) + 
  labs(title = "ROC curves for N0 network", x = "FPR", y = "TPR") +
  geom_line(size = 0.5) +
  scale_color_discrete(name = "Centrality",
                       breaks = levels(resultDF$Centrality),
                       labels = paste(c("DC", "BC","CC","EC","SGC"),", AUC: ",format(resultAUC,digits = 3),sep =''))
plot(myplot)

# ggsave(path = choose.dir(), filename = "N0_ROC.tiff", plot = myplot, device = "tiff", units = "in", width = 11, height = 8.5, dpi = 100)
