rm(list = ls())

library(muxViz)
library(igraph)
#> 
#> Attaching package: 'igraph'
#> The following objects are masked from 'package:stats':
#> 
#>     decompose, spectrum
#> The following object is masked from 'package:base':
#> 
#>     union
library(RColorBrewer)
library(viridis)
#> Loading required package: viridisLite
library(rgl)

source("~/Test_Git/muxViz/gui-old/muxVizGUI.R")


set.seed(1)
# Network setup
Layers <- 3
Nodes <- 200
layerCouplingStrength <- 1
networkOfLayersType <- "categorical"
isDirected <- F
layer.colors <- brewer.pal(8, "Set2")

pathInfomap <- "../src-exe/infomap-0.x/Infomap"


nodeTensor <- list()
g.list <- list()
plantedGroupsPerLayer <- 4
# matrix of the stochastic block model
block.matrix <- matrix(0.1 / Nodes, plantedGroupsPerLayer,
                       plantedGroupsPerLayer)
diag(block.matrix) <- 2 * log(Nodes) / Nodes
block.sizes <- rep(floor(Nodes / plantedGroupsPerLayer), plantedGroupsPerLayer)

for (l in 1:Layers) {
  #Generate the layers
  g.list[[l]] <- sample_sbm(Nodes, pref.matrix=block.matrix,
                            block.sizes=block.sizes)
  
  #Get the list of adjacency matrices which build the multiplex
  nodeTensor[[l]] <- get.adjacency(g.list[[l]])
}




lay <- layoutMultiplex(g.list, layout="fr", ggplot.format=F, box=T)

# Show the multiplex network
plot_multiplex3D(g.list, layer.layout=lay, layer.colors=layer.colors,
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels="auto", layer.labels.cex=1.5,
                 node.size.values="auto", node.size.scale=0.8,
                 show.aggregate=T)
# save the plot:
snapshot3d("multi_sbm.png", fmt = "png", width = 1024, height = 1024)
#> Warning in snapshot3d("../man/figures/multi_sbm.png", fmt = "png", width =
#> 1024, : webshot = TRUE requires the webshot2 package; using rgl.snapshot()
#> instead


