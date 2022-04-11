library(dplyr)
library(RCy3)
library(igraph)

PPI_base <- read.csv("~/Documents/List_genes_from_transcriptome/FIsInGene_020720_with_annotations.tsv", sep = "\t")

colnames(PPI_base) <- c("source", "target", "Annotation", "Direction", "Score")

deal_with_NA <- function(df, column, value){
  for (i in column) {
    data.table::set(df,which(is.na(df[[i]])), i,value)
  }
}

Make_specific_network <- function(Net, DEGs, msvip, logFC_threshold = 1.5, P.Value_DEG_threshold = 0.05, pval_msviper_threshold = 0.1){
  Genes <- dplyr::filter(DEGs, abs(logFC) > logFC_threshold, P.Value < P.Value_DEG_threshold) %>% .[,c(1,4,7)]
  TFs <- dplyr::filter(msvip$mrs_table, pval < pval_msviper_threshold) %>% .[,c(1,3,4)]
  Nodes_of_Interest <- c(Genes$ID, TFs$TF) %>% unique()
  Net <- dplyr::filter(Net, source %in% Nodes_of_Interest & target %in% Nodes_of_Interest)
  Net_features <- merge(Genes, TFs, by.x = "ID", by.y = "TF", all.x = T, all.y = T)
  Nodes_in_Net <- c(Net$source, Net$target) %>% .[. %ni%Net_features$ID] %>% unique()
  Nodes_in_Net <- data.frame(ID = Nodes_in_Net)
  Net_features <- merge(Net_features, Nodes_in_Net, by = "ID", all.x = T, all.y = T)
  Net_features <- dplyr::filter(Net_features, ID %in% c(Net$source, Net$target))
  deal_with_NA(Net_features, c("logFC", "nes"), 0)
  deal_with_NA(Net_features, c("P.Value", "pval"), 1)
  colnames(Net_features)[1] <- "id"
  igraph_net <- graph_from_data_frame(Net, vertices = Net_features, directed = T)
  eigen_centrality_result <- eigen_centrality(igraph_net, directed = F)$vector %>% as.data.frame()
  Net_features <- merge(Net_features, eigen_centrality_result, by.x = "id", by.y = 0)
  colnames(Net_features)[6] <- "Eigenvalue"
  list(Edges = Net, Nodes = Net_features)
}
