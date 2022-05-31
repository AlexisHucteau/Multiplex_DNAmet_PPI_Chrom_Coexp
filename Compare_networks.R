setwd("~/GitHub/Multiplex_DNAmet_PPI_Chrom_Coexp/")

source("~/Core_scripts/msviper_functions.R")
source("~/Core_scripts/core_functions.R")
source("~/GitHub/Multiplex_DNAmet_PPI_Chrom_Coexp/Make_koichi_factor.R")
# Aracne-AP networks

RNAseq <- read.csv("~/GitHub/Koichi_gene_expression_analyses_git/Koichi_gene_expression_analyses/DATA/RNAseq_parsed.csv", row.names = 1, header = T, check.names = F)
suppressPackageStartupMessages({
  library(dplyr)
  library(aracne.networks)
  library(GeneAnswers)
  library(org.Hs.eg.db)
  library(RCy3)
  library(stringr)
  library(viper)
  library(igraph)
  library(gridExtra)
  library(grid)
})

data(regulonlaml)
viper_regulons2dorothea <- function(r) {
  res <- r %>%
    purrr::map_df(
      .f = function(i) {
        tf_target <- i$tfmode %>%
          tibble::enframe(name = "target", value = "mor") %>%
          mutate(likelihood = i$likelihood)
      },
      .id = "tf"
    )
  return(res)
}

regulonaml <- viper_regulons2dorothea(regulonlaml)
regulonaml_SYMBOL <- data.frame("source" = GeneAnswers::getSymbols(regulonaml$tf, data = "org.Hs.eg"), 
                                "target" = GeneAnswers::getSymbols(regulonaml$target, data = "org.Hs.eg"),
                                "mor" = regulonaml$mor,
                                "likelihood" = regulonaml$likelihood)

data(dorothea_hs, package = "dorothea")
regulons_C = dorothea_hs %>%
  filter(confidence %in% c("A", "B", "C"))

regulons_C <- regulons_C[,c(1,3,2,4)]

regulons_B = dorothea_hs %>%
  filter(confidence %in% c("A", "B"))
         
regulons_B <- regulons_B[,c(1,3,2,4)]

allregulons = dorothea_hs

allregulons <- allregulons[,c(1,3,2,4)]

         
Bad_responder_aracne_network <- read.csv("output_Bad_responder/network.txt", sep = "\t") 
Good_responder_aracne_network <- read.csv("output_Good_responder/network.txt", sep = "\t")
Baseline_aracne_network <- read.csv("output_Baseline/network.txt", sep = "\t")
Good_n_Bad_aracne_network <- read.csv("output_Good_n_Bad_responder/network.txt", sep = "\t")
Relapse_n_Good_aracn_network <- read.csv("output_Relapse_n_Good/network.txt", sep = "\t")
All_sample_aracn_network <- read.csv("output/network.txt", sep = "\t")

ARACNe_networks <- list(Bad_responder_aracne_network, Baseline_aracne_network, Good_n_Bad_aracne_network, Relapse_n_Good_aracn_network, All_sample_aracn_network, regulonaml_SYMBOL, regulons_B, regulons_C, allregulons)
names(ARACNe_networks) <- c("Bad_responder_aracne_network", "Baseline_aracne_network", "Good_n_Bad_aracne_network", "Relapse_n_Good_aracn_network", "All_sample_aracn_network", "regulonaml_SYMBOL", "regulons_B", "regulons_C", "allregulons")

ARACNe_networks_cytos <- lapply(names(ARACNe_networks), function(net){
  colnames(ARACNe_networks[[net]]) <- c("source", "target", "MI", "pvalue")
  createNetworkFromDataFrames(edges = ARACNe_networks[[net]], title = net, collection = "ARACNe network6")
})

NR_R_ARACNe_networks <- c("Bad_responder_aracne_network", "Baseline_aracne_network", "Good_n_Bad_aracne_network", "All_sample_aracn_network", "regulonaml_SYMBOL", "regulons_B", "regulons_C", "allregulons")

ref_R_B <- Factor_R_OR_NR_B == "R.B"
ref_NR_B <- Factor_R_OR_NR_B == "NR.B"

ref_REL_POST <- Factor_R_OR_NR_B == "OR.REL" | Factor_R_OR_NR_B == "R.REL"

source("network_functions.R")
RNAseq_diff_exp <- readRDS("~/GitHub/Koichi_gene_expression_analyses_git/Koichi_gene_expression_analyses/Results/Tables/RNAseq_diff_gene_expression_analysis.rds")

NR_R_msviper_networks <- lapply(NR_R_ARACNe_networks, function(net){
  colnames(ARACNe_networks[[net]]) <- c("tf", "target", "MI", "pvalue")
  ms <- run_msviper(RNAseq, ARACNe_networks[[net]], use_aracne = T, ref_R_B, ref_NR_B,  "R", "NR", minsize = 4, ges.filter=T)
  colnames(ms$regulons) <- c("source", "target", "mor", "likelihood", "state")
  # msviper_network <- Make_specific_network(Net = ms$regulons, DEGs = RNAseq_diff_exp$R_OR_NR_B$`NR.B-R.B`, msvip = ms, title = net, collection = "NR_vs_R", P_P = F, PPI = F)
  list(ms)
})

names(NR_R_msviper_networks) <- NR_R_ARACNe_networks

REL_R_ARACNe_networks <- c("Baseline_aracne_network", "All_sample_aracn_network", "regulonaml_SYMBOL", "regulons_B", "regulons_C", "allregulons", "Relapse_n_Good_aracn_network")

REL_R_msviper_networks <- lapply(REL_R_ARACNe_networks, function(net){
  colnames(ARACNe_networks[[net]]) <- c("tf", "target", "MI", "pvalue")
  ms <- run_msviper(RNAseq, ARACNe_networks[[net]], use_aracne = T, ref_R_B, ref_REL_POST,  "R", "REL", minsize = 4, ges.filter=T)
  colnames(ms$regulons) <- c("source", "target", "mor", "likelihood", "state")
  msviper_network <- Make_specific_network(Net = ms$regulons, DEGs = RNAseq_diff_exp$R_OR_NR_B$`R.REL-R.B`, msvip = ms, title = net, collection = "REL_vs_R", P_P = F, PPI = F)
  list(ms,
       msviper_network)
})

names(REL_R_msviper_networks) <- REL_R_ARACNe_networks


Degree_networks <- lapply(names(ARACNe_networks), function(net){
  igraphe_net <- graph_from_data_frame(ARACNe_networks[[net]])
  deg <- igraph::degree(igraphe_net)
  list("igraph" = igraphe_net, "degree" = deg)
})

names(Degree_networks) <- c("Bad_responder_aracne_network", "Baseline_aracne_network", "Good_n_Bad_aracne_network", "Relapse_n_Good_aracn_network", "All_sample_aracn_network", "regulonaml_SYMBOL", "regulons_B", "regulons_C", "allregulons")

enrichment_degree <- lapply(NR_R_ARACNe_networks, function(net){
  merge(NR_R_msviper_networks[[net]]$ms, )
})

viper_plot <- lapply(names(NR_R_msviper_networks), function(ARACNe_net){
  png(file = paste0("~/GitHub/Multiplex_DNAmet_PPI_Chrom_Coexp/ms_vipers_", ARACNe_net, "_NR_R.png"), width = 960, height = 540)
  plot(NR_R_msviper_networks[[ARACNe_net]][[1]]$mrs)
  dev.off()
})

names(viper_plot) <- NR_R_ARACNe_networks

png(file = paste0("~/GitHub/Multiplex_DNAmet_PPI_Chrom_Coexp/ms_vipers_NR_R.png"), width = 1920, height = 1080)
grid.arrange(plot(NR_R_msviper_networks[["Bad_responder_aracne_network"]][[1]]$mrs), plot(NR_R_msviper_networks[["Baseline_aracne_network"]][[1]]$mrs), plot(NR_R_msviper_networks[["Good_n_Bad_aracne_network"]][[1]]$mrs), plot(NR_R_msviper_networks[["All_sample_aracn_network"]][[1]]$mrs),
             ncol=3,
             top = textGrob('Chas based on response',
                            just = c('center'),
                            gp = gpar(fontsize = 32)))
dev.off()

png(file = paste0("~/GitHub/Multiplex_DNAmet_PPI_Chrom_Coexp/ms_vipers_REL_R_GoodnRelapse.png"), width = 1920, height = 1080)
plot(REL_R_msviper_networks$Relapse_n_Good_aracn_network[[1]]$mrs)
dev.off()
png(file = paste0("~/GitHub/Multiplex_DNAmet_PPI_Chrom_Coexp/ms_vipers_REL_R_allsamples.png"), width = 1920, height = 1080)
plot(REL_R_msviper_networks$All_sample_aracn_network[[1]]$mrs)
dev.off()
png(file = paste0("~/GitHub/Multiplex_DNAmet_PPI_Chrom_Coexp/ms_vipers_REL_R_regulons_B.png"), width = 1920, height = 1080)
plot(REL_R_msviper_networks$regulons_B[[1]]$mrs)
dev.off()


viper_plot <- lapply(names(REL_R_msviper_networks), function(ARACNe_net){
  png(file = paste0("~/GitHub/Multiplex_DNAmet_PPI_Chrom_Coexp/ms_vipers_", ARACNe_net, "_REL_R.png"), width = 960, height = 540)
  plot(REL_R_msviper_networks[[ARACNe_net]][[1]]$mrs)
  dev.off()
})