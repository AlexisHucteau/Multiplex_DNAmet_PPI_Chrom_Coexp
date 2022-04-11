library(dplyr)
library(RCy3)

prepare_pchic <- function(cell_lines = "all", minimum_interaction = 5){
  load("DATA/pchic.RData")
  if (length(cell_lines) >= 1){
    cell_lines = c("Mon", "Mac0", "Mac1", "Mac2", "Neu", "MK", "EP", "Ery", "FoeT", "nCD4", "tCD4", "aCD4", "naCD4", "nCD8", "tCD8", "nB", "tB")
  }
  pchic <- data.frame(pchic[rowSums(pchic[,cell_lines] >= minimum_interaction) >= 1, 1:10]) %>% na.omit(.)
  colnames(pchic)[c(1:5, 6:10)] <- rep(c("chr", "start", "end", "ID", "Name"), 2)
  return(pchic)
}

pchic <- prepare_pchic(cell_lines = c("Mon", "Mac1", "Mac0", "Mac2", "MK", "Ery", "EP"))
pchic_bed <- unique(rbind(pchic[, c(1:4, 5)], pchic[, c(6:9, 10)]))

Pchic_edges <- pchic[,c(4,9)]
colnames(Pchic_edges) <- c("source", "target")
Pchic_nodes <- pchic_bed[,c(4,5)]
Pchic_nodes$Promoter <- Pchic_nodes$ID %in% Pchic_edges$source
colnames(Pchic_nodes) <- c("id", "Gene_name", "Is_Promoter")

