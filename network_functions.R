deal_with_NA <- function(df, column, value){
  for (i in column) {
    data.table::set(df,which(is.na(df[[i]])), i,value)
  }
}

Make_specific_network <- function(Net, DEGs, msvip, title, collection, P_P = T, PPI=T, logFC_threshold = 1.5, P.Value_DEG_threshold = 0.05, pval_msviper_threshold = 0.1){
  Genes <- dplyr::filter(DEGs, abs(logFC) > logFC_threshold, P.Value < P.Value_DEG_threshold) %>% .[,c(1,4,7)]
  TFs <- dplyr::filter(msvip$mrs_table, pval < pval_msviper_threshold) %>% .[,c(1,3,4)]
  Nodes_of_Interest <- c(Genes$ID, TFs$TF) %>% unique()
  if (P_P){
    Net <- dplyr::filter(Net, source %in% Nodes_of_Interest & target %in% Nodes_of_Interest)
  }else{
    Net <- dplyr::filter(Net, source %in% Nodes_of_Interest | target %in% Nodes_of_Interest)
  }
  Net_features <- merge(DEGs, msvip$mrs_table, by.x = "ID", by.y = "TF", all.x = T, all.y = T)
  Nodes_in_Net <- c(Net$source, Net$target) %>% .[. %ni%Net_features$ID] %>% unique()
  Nodes_in_Net <- data.frame(ID = Nodes_in_Net)
  Net_features <- merge(Net_features, Nodes_in_Net, by = "ID", all.x = T, all.y = T)
  Net_features <- dplyr::filter(Net_features, ID %in% c(Net$source, Net$target))
  deal_with_NA(Net_features, c("logFC", "nes"), 0)
  deal_with_NA(Net_features, c("P.Value", "pval"), 1)
  colnames(Net_features)[1] <- "id"
  igraph_net <- graph_from_data_frame(Net, vertices = Net_features, directed = T)
  eigen_centrality_result <- eigen_centrality(igraph_net, directed = F)$vector %>% as.data.frame()
  Net_features <- Net_features[,c("id", "logFC", "P.Value", "nes", "pval")]
  Net_features <- merge(Net_features, eigen_centrality_result, by.x = "id", by.y = 0)
  colnames(Net_features)[6] <- "Eigen_centrality"
  Net_features$Type <- ifelse(Net_features$nes == 0, "Gene", "TF")
  createNetworkFromDataFrames(Net_features, Net, title, collection)
  if(PPI){
    Modify_Cytoscape_network_PPI(Net, Net_features, title, collection)
  }else{
    Modify_Cytoscape_network_Aracn(Net, Net_features, title, collection)
  }
  list(Edges = Net, Nodes = Net_features)
  
}


Modify_Cytoscape_network_PPI <- function(e, n, title, collection){
  defaults <- list(NODE_SHAPE="diamond",
                   NODE_SIZE=10,
                   EDGE_TRANSPARENCY=120,
                   NODE_LABEL_POSITION="c,c,c,0.00,0.00")
  nodeLabels <- mapVisualProperty(visual.prop = "node label", table.column = 'name', mapping.type = 'p')
  createVisualStyle(title, defaults, list(nodeLabels))
  setVisualStyle(title)
  setNodeShapeMapping(table.column = "Type", 
                      table.column.values = c("TF", "Gene"),
                      shapes = c('ELLIPSE', 'RECTANGLE'),
                      style.name = title)
  
  setNodeColorMapping(table.column = 'logFC',
                      table.column.values = c(min(n$logFC), 0.0, max(n$logFC)),
                      colors = c('#0000FF', '#FFFFFF', '#FF0000'),
                      style.name = title)
  
  setNodeFillOpacityMapping(table.column = 'P.Value', 
                            table.column.values = c(0, 0.1, 1), 
                            opacities = c(255, 200, 100), 
                            style.name = title)
  
  setNodeSizeMapping (table.column = 'Eigen_centrality', 
                      table.column.values = c(0, 1), 
                      sizes = c(30, 200), 
                      style.name = title)
  
  setNodeFontSizeMapping(table.column = 'Eigen_centrality', 
                         table.column.values = c(0, 1), 
                         sizes = c(15, 75), 
                         style.name = title)
  
  setEdgeTargetArrowShapeMapping(table.column = 'Direction',
                                 table.column.values = names(table(e$Direction)),
                                 shapes = c('NONE', 'ARROW', 'T', 'NONE', 'ARROW', 'T', 'NONE', 'ARROW'), 
                                 style.name = title)
  
  createDegreeFilter(filter.name = "single", criterion = c(0,0))
  deleteSelectedNodes()
  
  createColumnFilter(filter.name = "NES_DOWN", column = "nes", criterion = 0, type = "nodes", predicate = "LESS_THAN")
  if (!is.na(getSelectedNodes())){
    setNodeColorBypass(getSelectedNodes(), new.colors = "#C411FF")
    clearSelection()
  }
  createColumnFilter(filter.name = "NES_UP", column = "nes", criterion = 0, type = "nodes", predicate = "GREATER_THAN")
  if (!is.na(getSelectedNodes())){
    setNodeColorBypass(getSelectedNodes(), new.colors = "#5AFF00")  
    clearSelection()
  }
  createColumnFilter(filter.name = "TF_activity_pval", column = "pval", criterion = c(0.1,0.99), type = "nodes", predicate = "BETWEEN")
  if (!is.na(getSelectedNodes())){
    setNodeColorBypass(getSelectedNodes(), new.colors = "#BBBBBB")
    clearSelection()
  }
  layoutNetwork()
  exportImage(paste0(collection, "_", title), 'SVG', zoom=200)
}


Modify_Cytoscape_network_Aracn <- function(e, n, title, collection){
  defaults <- list(NODE_SHAPE="diamond",
                   NODE_SIZE=10,
                   EDGE_TRANSPARENCY=120,
                   NODE_LABEL_POSITION="c,c,c,0.00,0.00")
  nodeLabels <- mapVisualProperty(visual.prop = "node label", table.column = 'name', mapping.type = 'p')
  createVisualStyle(title, defaults, list(nodeLabels))
  setVisualStyle(title)
  setNodeShapeMapping(table.column = "Type", 
                      table.column.values = c("TF", "Gene"),
                      shapes = c('ELLIPSE', 'RECTANGLE'),
                      style.name = title)
  
  setNodeColorMapping(table.column = 'logFC',
                      table.column.values = c(min(n$logFC), 0.0, max(n$logFC)),
                      colors = c('#0000FF', '#FFFFFF', '#FF0000'),
                      style.name = title)
  
  setNodeFillOpacityMapping(table.column = 'P.Value', 
                            table.column.values = c(0, 0.1, 1), 
                            opacities = c(255, 200, 100), 
                            style.name = title)
  
  setNodeSizeMapping (table.column = 'Eigen_centrality', 
                      table.column.values = c(0, 1), 
                      sizes = c(30, 200), 
                      style.name = title)
  
  setNodeFontSizeMapping(table.column = 'Eigen_centrality', 
                         table.column.values = c(0, 1), 
                         sizes = c(15, 75), 
                         style.name = title)
  
  setEdgeLineWidthMapping(table.column = 'mor',
                          table.column.values = c(min(e$mor), 0.0, max(e$mor)),
                          widths = c(10, 2, 10), 
                          style.name = title)
  
  setEdgeColorMapping(table.column = 'mor',
                      table.column.values = c(min(e$mor), 0.0, max(e$mor)),
                      colors = c('#0000FF', "#555555", '#FF0000'), 
                      style.name = title)
  
  setEdgeOpacityMapping(table.column = 'mor',
                        table.column.values = c(min(e$mor), 0.0, max(e$mor)),
                        opacities = c(255, 100, 255), 
                        style.name = title)
  
  createDegreeFilter(filter.name = "single", criterion = c(0,0))
  deleteSelectedNodes()
  
  createColumnFilter(filter.name = "NES_DOWN", column = "nes", criterion = 0, type = "nodes", predicate = "LESS_THAN")
  setNodeColorBypass(getSelectedNodes(), new.colors = "#C411FF")
  createColumnFilter(filter.name = "NES_UP", column = "nes", criterion = 0, type = "nodes", predicate = "GREATER_THAN")
  setNodeColorBypass(getSelectedNodes(), new.colors = "#5AFF00")  
  createColumnFilter(filter.name = "TF_activity_pval", column = "pval", criterion = c(0.1,0.99), type = "nodes", predicate = "BETWEEN")
  setNodeColorBypass(getSelectedNodes(), new.colors = "#BBBBBB")
  clearSelection()
  layoutNetwork()
  exportImage(paste0(collection, "_", title), 'SVG', zoom=200)
}

Make_pchic_genes_specific <- function(Pchic_net = pchic_prom, DEGs, TFs, title, collection,  P_P = T, logFC_threshold = 1.5, P.Value_DEG_threshold = 0.05, pval_msviper_threshold = 0.1){
  Genes <- dplyr::filter(DEGs, abs(logFC) > logFC_threshold, P.Value < P.Value_DEG_threshold) %>% .[,c(1,4,7)]
  TFs <- dplyr::filter(TFs$mrs_table, pval < pval_msviper_threshold) %>% .[,c(1,3,4)]
  Promoter_of_interest <- c(Genes$ID, TFs$TF) %>% unique()
  if (P_P){
    Pchic_net <- dplyr::filter(pchic_prom, Name1 %in% Promoter_of_interest & Name2 %in% Promoter_of_interest)
  }else{
    Pchic_net <- dplyr::filter(pchic_prom, Name1 %in% Promoter_of_interest | Name2 %in% Promoter_of_interest)
  }
  Pchic_net$Name2 <- ifelse(Pchic_net$Name2 == ".", Pchic_net$ID2, Pchic_net$Name2)
  Pchic_net <- unique(Pchic_net[,c(5,10)])
  Pchic_net_features <- c(Pchic_net$Name1, Pchic_net$Name2) %>% unique()
  Pchic_net_features <- data.frame(ID = Pchic_net_features, na = rep(0, length(Pchic_net_features)))
  Pchic_net_features <- merge(Pchic_net_features, Genes, by.x = "ID", by.y = "ID", all.x = T)
  Pchic_net_features <- merge(Pchic_net_features, TFs, by.x = "ID", by.y = "TF", all.x = T)
  Pchic_net_features <- Pchic_net_features[,c(1,3:6)]
  deal_with_NA(Pchic_net_features, c("logFC", "nes"), 0)
  deal_with_NA(Pchic_net_features, c("P.Value", "pval"), 1)
  colnames(Pchic_net_features)[1] <- "id"
  igraph_net <- graph_from_data_frame(Pchic_net, vertices = Pchic_net_features, directed = T)
  eigen_centrality_result <- eigen_centrality(igraph_net, directed = F)$vector %>% as.data.frame()
  Pchic_net_features <- merge(Pchic_net_features, eigen_centrality_result, by.x = "id", by.y = 0)
  colnames(Pchic_net_features)[6] <- "Eigen_centrality"
  colnames(Pchic_net) <- c("source", "target")
  Pchic_net_features$Type <- ifelse(Pchic_net_features$nes == 0, "Gene", "TF")
  createNetworkFromDataFrames(Pchic_net_features, Pchic_net, title, collection)
  Modify_Cytoscape_network_Pchic_DEGs(Pchic_net, Pchic_net_features, title, collection)
  list(Edges = Pchic_net, Nodes = Pchic_net_features)
}

Modify_Cytoscape_network_Pchic_DEGs <- function(e, n, title, collection){
  defaults <- list(NODE_SHAPE="diamond",
                   NODE_SIZE=10,
                   EDGE_TRANSPARENCY=120,
                   NODE_LABEL_POSITION="c,c,c,0.00,0.00")
  nodeLabels <- mapVisualProperty(visual.prop = "node label", table.column = 'name', mapping.type = 'p')
  createVisualStyle(title, defaults, list(nodeLabels))
  setVisualStyle(title)
  setNodeShapeMapping(table.column = "Type", 
                      table.column.values = c("TF", "Gene"),
                      shapes = c('ELLIPSE', 'RECTANGLE'),
                      style.name = title)
  
  setNodeColorMapping(table.column = 'logFC',
                      table.column.values = c(min(n$logFC), 0.0, max(n$logFC)),
                      colors = c('#0000FF', '#FFFFFF', '#FF0000'),
                      style.name = title)
  
  setNodeFillOpacityMapping(table.column = 'P.Value', 
                            table.column.values = c(0, 0.1, 1), 
                            opacities = c(255, 200, 100), 
                            style.name = title)
  
  createDegreeFilter(filter.name = "single", criterion = c(0,0))
  deleteSelectedNodes()
  
  createColumnFilter(filter.name = "NES_DOWN", column = "nes", criterion = 0, type = "nodes", predicate = "LESS_THAN")
  if (!is.na(getSelectedNodes())){
    setNodeColorBypass(getSelectedNodes(), new.colors = "#C411FF")
  }
  createColumnFilter(filter.name = "NES_UP", column = "nes", criterion = 0, type = "nodes", predicate = "GREATER_THAN")
  if (!is.na(getSelectedNodes())){
    setNodeColorBypass(getSelectedNodes(), new.colors = "#5AFF00")  
  }
  createColumnFilter(filter.name = "TF_activity_pval", column = "pval", criterion = c(0.1,0.99), type = "nodes", predicate = "BETWEEN")
  if (!is.na(getSelectedNodes())){
    setNodeColorBypass(getSelectedNodes(), new.colors = "#BBBBBB")
  }
  clearSelection()
  layoutNetwork()
  exportImage(paste0(collection, "_", title), 'SVG', zoom=200)
}

Make_pchic_methylation_network <- function(DMPs, DMRs, pchic_net = pchic, P_P = T, title, collection){
  colnames(pchic_net) <- c("chr1", "start1", "end1", "ID1", "Name1", "chr2", "start2", "end2", "ID2", "Name2")
  DMPs <- na.omit(DMPs)
  DMRs <- na.omit(DMRs)
  DMPs <- DMPs[,c(1,2,5,23)]
  DMRs <- DMRs[,c(1, 7, 14, 19)]
  Methylation_fragments <- c(DMPs$ID, DMRs$ID) %>% unique()
  if (P_P){
    pchic_net <- dplyr::filter(pchic_net, ID1 %in% Methylation_fragments & ID2 %in% Methylation_fragments)
  }else{
    pchic_net <- dplyr::filter(pchic_net, ID1 %in% Methylation_fragments | ID2 %in% Methylation_fragments)
  }
  pchic_net <- pchic_net[,c(4,5,9,10)]
  pchic_net_final <- pchic_net[,c(1,3)]
  colnames(pchic_net) <- c("ID", "Name", "ID", "Name")
  pchic_net_features <- rbind(pchic_net[,c(1,2)], pchic_net[,c(3,4)])
  pchic_features <- merge(pchic_net_features, DMPs, by.x = "ID", by.y = "ID", all.x = T)
  pchic_features <- merge(pchic_features, DMRs, by.x = "ID", by.y = "ID", all.x = T)
  deal_with_NA(pchic_features, c("rowname.x", "rowname.y"), "")
  deal_with_NA(pchic_features, c("logFC", "value"), 0)
  deal_with_NA(pchic_features, c("p.value", "P.Value"), 1)
  colnames(pchic_features) <- c("id", "Name", "cpg", "logFC_cpg", "P.Value_cpg", "DMR", "logFC_DMR", "P.Value_DMR")
  pchic_features <- unique(pchic_features)
  pchic_features <- pchic_features[!duplicated(pchic_features$id),]
  igraph_net <- graph_from_data_frame(pchic_net_final, vertices = pchic_features, directed = T)
  eigen_centrality_result <- eigen_centrality(igraph_net, directed = F)$vector %>% as.data.frame()
  Pchic_net_features <- merge(pchic_features, eigen_centrality_result, by.x = "id", by.y = 0)
  colnames(Pchic_net_features)[9] <- "Eigenvalue"
  colnames(pchic_net_final) <- c("source", "target")
  Pchic_net_features$id <- as.character(Pchic_net_features$id)
  Pchic_net_features$Fragment_type <- paste0(Pchic_net_features$cpg, Pchic_net_features$DMR, Pchic_net_features$Name)
  Pchic_net_features$Fragment_type <- sapply(Pchic_net_features$Fragment_type, function(frag){
    if(!stringr::str_detect(frag, "[:punct:]$")){
      "Gene"
    }else if(stringr::str_detect(frag, "DMR")){
      "DMR"
    }else if(stringr::str_detect(frag, "cg")){
      "cpg"
    }else{
      "fragment"
    }
  })
  pchic_net_final$source <- as.character(pchic_net_final$source)
  pchic_net_final$target <- as.character(pchic_net_final$target)
  createNetworkFromDataFrames(Pchic_net_features, pchic_net_final, title, collection)
  Modify_Cytoscape_network_Pchic_Methyls(Pchic_net, Pchic_net_features, title, collection)
  
  
  list(Edges = pchic_net_final, Nodes = Pchic_net_features)
}

Modify_Cytoscape_network_Pchic_Methyls <- function(e, n, title, collection){
  defaults <- list(NODE_SHAPE="diamond",
                   NODE_SIZE=10,
                   EDGE_TRANSPARENCY=120,
                   NODE_LABEL_POSITION="c,c,c,0.00,0.00")
  nodeLabels <- mapVisualProperty(visual.prop = "node label", table.column = 'name', mapping.type = 'p')
  createVisualStyle(title, defaults, list(nodeLabels))
  setVisualStyle(title)
  setNodeShapeMapping(table.column = "Fragment_type", 
                      table.column.values = c("DMR", "cpg", "Gene", "fragment"),
                      shapes = c('RECTANGLE', 'VEE', 'ELLIPSE', 'DIAMOND'),
                      style.name = title)
  
  setNodeColorMapping(table.column = 'logFC_cpg',
                      table.column.values = c(min(n$logFC_cpg), 0.0, max(n$logFC_cpg)),
                      colors = c('#0000FF', '#FFFFFF', '#FF0000'),
                      style.name = title)
  
  setNodeFillOpacityMapping(table.column = 'P.Value_cpg', 
                            table.column.values = c(0, 0.1, 1), 
                            opacities = c(255, 200, 100), 
                            style.name = title)
  
  setNodeSizeMapping(table.column = 'Fragment_type', 
                     table.column.values = c("DMR", "cpg", "Gene", "fragment"), 
                     sizes = c(40, 20, 60, 1), 
                     style.name = title, mapping.type = 'd')
  
  setNodeFontSizeMapping(table.column = 'Fragment_type', 
                         table.column.values = c("DMR", "cpg", "Gene", "fragment"), 
                         sizes = c(40, 20, 60, 1), 
                         style.name = title, mapping.type = 'd')
  
  createDegreeFilter(filter.name = "single", criterion = c(0,0))
  deleteSelectedNodes()
  
  createColumnFilter(filter.name = "Gene", column = "Fragment_type", criterion = "Gene", type = "nodes", predicate = "IS")
  if (!is.na(getSelectedNodes())){
    setNodeSizeBypass(getSelectedNodes(), new.sizes = 60)
    setNodeFontSizeBypass(getSelectedNodes(), new.sizes = 60)
  }
  clearSelection()
  layoutNetwork()
  exportImage(paste0(collection, "_", title), 'SVG', zoom=200)
}