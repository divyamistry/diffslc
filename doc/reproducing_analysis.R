getwd()
setwd('~/Desktop/reproducible/')

source('src/pre-prep/get_Tu_etal_data.R')
source('src/pre-prep/get_Guan_etal_data.R')
source('src/pre-prep/process_DIP_20150101_data.R')
source('src/pre-prep/DIP_to_Affy_Mapping.R')
source('src/pre-prep/make_network.R')
source('src/graph_statistics_gse3431.R')
source('src/graph_statistics_gse3076.R')

require(igraph)
list.graph.attributes(graph = hybridG.gse3431) #for DIP + Tu et.al.
list.vertex.attributes(graph = hybridG.gse3431)
list.edge.attributes(graph = hybridG.gse3431)

# Plots have their own customization so they're not treated nicely
# like having their own functions. Perhaps a generalized function
# for this might help, but I think the linear code in somePlots.R
# is much easier to read and follow along the way it is now.
hybridG <- hybridG.gse3431
graphRankingStats <- graphRankingStats.gse3431
source('src/somePlots_N0.R') # N0 doesnt use gene expression, so this can be done using any network.
source('src/somePlots_NT.R')

hybridG <- hybridG.gse3076
graphRankingStats <- graphRankingStats.gse3076
source('src/somePlots_NF.R')

source('src/somePlots_NT-NF-Comparison.R') # data for this is already processed and available.

rm(hybridG, graphRankingStats)

source('src/getOrderedCummulativeCounts.R')
(df <- data.frame(Col1 = rnorm(20), Col2 = runif(20)/rnorm(20,mean = 20), Response = (sample(2000,20) > 1000)))
getOrderedCummulativeCounts(df = df, countColumn = "Response", decreasing = T)

# a test run on my Early-2011 MBP took ~105 minutes.

