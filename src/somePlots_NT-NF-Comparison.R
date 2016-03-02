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

#' PR curve for degree vs NT3 and NF3.
plot(
  performance(
    prediction.obj = prediction(predictions = graphRankingStats$degree, labels = graphRankingStats$essential), 
    measure = "prec", x.measure = "rec"), 
  main = "Comparison between DC and DiffSLc", 
  add = F, col = "black"
)

plot(
  performance(
    prediction.obj = prediction(predictions = graphRankingStats$eigcent, labels = graphRankingStats$essential), 
    measure = "prec", x.measure = "rec"), 
  add = T, col = "#919191"
)

plot(
  performance(
    prediction.obj = prediction(predictions = graphRankingStats$NT3_diffslc, labels = graphRankingStats$essential), 
    measure = "prec", x.measure = "rec"), 
  add = T, col = "orange"
)

plot(
  performance(
    prediction.obj = prediction(predictions = graphRankingStats$NF3_diffslc, labels = graphRankingStats$essential), 
    measure = "prec", x.measure = "rec"), 
  add = T, col = "blue"
)

legend(x = "topright", legend = c("Eigen vector (EC)", "Degree (DC)", "DiffSLc NT3", "DiffSLc NF3"), 
       col = c("#919191", "black","orange","blue"), lty = c(1,1,1,1), cex = 1)

######## 
# plot(sort(graphRankingStats$eigcent), type="l", xlab = "", ylab = "ev centrality")
# plot(sort(graphRankingStats$betweenness), type="l", xlab = "", ylab = "betweenness")
# plot(sort(graphRankingStats$mdc2Max7), type="l", xlab = "", ylab = "mdc2Max7")