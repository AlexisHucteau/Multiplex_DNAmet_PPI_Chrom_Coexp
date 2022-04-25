library(dplyr)
setwd("GitHub/Multiplex_DNAmet_PPI_Chrom_Coexp/")

# Aracne-AP networks

RNAseq <- read.csv("~/GitHub/Koichi_gene_expression_analyses_git/Koichi_gene_expression_analyses/DATA/RNAseq_parsed.csv", row.names = 1, header = T, check.names = F)

regulonaml <- viper_regulons2dorothea(regulonlaml)
regulonaml_SYMBOL <- data.frame("source" = GeneAnswers::getSymbols(regulonaml$tf, data = "org.Hs.eg"), 
                                "target" = GeneAnswers::getSymbols(regulonaml$target, data = "org.Hs.eg"),
                                "mor" = regulonaml$mor,
                                "likelihood" = regulonaml$likelihood)

data(dorothea_hs, package = "dorothea")
regulons_C = dorothea_hs %>%
  filter(confidence %in% c("A", "B", "C"))

regulons_B = dorothea_hs %>%
  filter(confidence %in% c("A", "B"))
         
allregulons = dorothea_hs
         
Bad_responder_aracne_network <- read.csv("output_Bad_responder/network.txt", sep = "\t") 
Good_responder_aracne_network <- read.csv("output_Good_responder/network.txt", sep = "\t")
Baseline_aracne_network <- read.csv("output_Baseline/network.txt", sep = "\t")
Good_n_Bad_aracne_network <- read.csv("output_Good_n_Bad_responder/network.txt", sep = "\t")
Relapse_n_Good_aracn_network <- read.csv("output_Relapse_n_Good/network.txt", sep = "\t")
All_sample_aracn_network <- read.csv("output/network.txt", sep = "\t")

ARACNe_networks <- list(Bad_responder_aracne_network, Baseline_aracne_network, Good_n_Bad_aracne_network, Relapse_n_Good_aracn_network, All_sample_aracn_network, regulonaml_SYMBOL, regulons_B, regulons_C, allregulons)
names(ARACNe_networks) <- c("Bad_responder_aracne_network", "Baseline_aracne_network", "Good_n_Bad_aracne_network", "Relapse_n_Good_aracn_network", "All_sample_aracn_network", "regulonaml_SYMBOL", "regulons_B", "regulons_B", "allregulons")

ARACNe_networks_cytos <- lapply(names(ARACNe_networks), function(net){
  colnames(ARACNe_networks[[net]]) <- c("source", "target", "MI", "pvalue")
  createNetworkFromDataFrames(edges = ARACNe_networks[[net]], title = net, collection = "ARACNe network")
})

NR_R_ARACNe_networks <- c("Bad_responder_aracne_network", "Baseline_aracne_network", "Good_n_Bad_aracne_network", "All_sample_aracn_network", "regulonaml_SYMBOL", "regulons_B", "regulons_B", "allregulons")

ref_R_B <- Factor_R_OR_NR_B == "R.B"
ref_NR_B <- Factor_R_OR_NR_B == "NR.B"

ref_REL_POST <- Factor_R_OR_NR_B == "OR.REL" | Factor_R_OR_NR_B == "R.REL"


NR_R_msviper_networks <- lapply(NR_R_ARACNe_networks, function(net){
  colnames(ARACNe_networks[[net]]) <- c("tf", "target", "MI", "pvalue")
  ms <- run_msviper(RNAseq, ARACNe_networks[[net]], use_aracne = T, ref_R_B, ref_NR_B,  "R", "NR", minsize = 4, ges.filter=T)
  colnames(ms$regulons) <- c("source", "target", "mor", "likelihood", "state")
  msviper_network <- Make_specific_network(Net = ms$regulons, DEGs = RNAseq_diff_gene_expression_analysis$R_OR_NR_B$`NR.B-R.B`, msvip = ms, title = net, collection = "NR_vs_R", P_P = F, PPI = F)
  list(ms,
       msviper_network)
})



colnames(Bad_responder_aracne_network)[c(1,2)] <- c("source", "target")
colnames(Good_responder_aracne_network)[c(1,2)] <- c("source", "target")
colnames(Baseline_aracne_network)[c(1,2)] <- c("source", "target")
colnames(Good_n_Bad_aracne_network)[c(1,2)] <- c("source", "target")

# msviper aracn networks

NR_R_msviper <- readRDS("~/GitHub/Koichi_gene_expression_analyses_git/Koichi_gene_expression_analyses/Results/Tables/NR_R_msviper.rds")
RNAseq_diff_exp <- readRDS("~/GitHub/Koichi_gene_expression_analyses_git/Koichi_gene_expression_analyses/Results/Tables/RNAseq_diff_gene_expression_analysis.rds")

source("~/Core_scripts/msviper_functions.R")
suppressPackageStartupMessages({
  library(viper)
  library(aracne.networks)
  library(GeneAnswers)
  library(org.Hs.eg.db)
  library(RCy3)
  library(igraph)
  library(stringr)
})

data(regulonlaml)

Bad_responder_aracne_igraph <- igraph::graph_from_data_frame(Bad_responder_aracne_network)

Bad_responder_aracne_degree <- as.data.frame(degree(Bad_responder_aracne_igraph))

Baseline_aracne_igraph <- igraph::graph_from_data_frame(Baseline_aracne_network)

Baseline_aracne_degree <- as.data.frame(degree(Baseline_aracne_igraph))

Good_n_Bad_aracne_igraph <- igraph::graph_from_data_frame(Good_n_Bad_aracne_network)

Good_n_Bad_aracne_degree <- as.data.frame(degree(Good_n_Bad_aracne_igraph))


ref_R_B <- Factor_R_OR_NR_B == "R.B"
ref_NR_B <- Factor_R_OR_NR_B == "NR.B"

ref_REL_POST <- Factor_R_OR_NR_B == "OR.REL" | Factor_R_OR_NR_B == "R.REL"

regulon_Bad_responder <- Bad_responder_aracne_network
colnames(regulon_Bad_responder)[1:2] <- c("tf", "target")

NR_R_Bad_aracn_msviper <- run_msviper(RNAseq, regulon_Bad_responder, use_aracne = T, ref_R_B, ref_NR_B,  "R", "NR", minsize = 4, ges.filter=T)
colnames(NR_R_Bad_aracn_msviper$regulons)[1:2] <- c("source", "target")

NR_R_Bad_aracn_msviper_network <- Make_specific_network(Net = NR_R_Bad_aracn_msviper$regulons, DEGs = RNAseq_diff_gene_expression_analysis$R_OR_NR_B$`NR.B-R.B`, msvip = NR_R_Bad_aracn_msviper, title = "Bad_aracn_network", collection = "NR_vs_R", P_P = F, PPI = F)

plot(NR_R_Bad_aracn_msviper$mrs, cex = .7)


regulon_Good_n_Bad_responder <- Good_n_Bad_aracne_network
colnames(regulon_Good_n_Bad_responder)[1:2] <- c("tf", "target")

NR_R_Good_n_Bad_aracn_msviper <- run_msviper(RNAseq, regulon_Good_n_Bad_responder, use_aracne = T, ref_R_B, ref_NR_B,  "R", "NR", minsize = 4, ges.filter=T)
colnames(NR_R_Good_n_Bad_aracn_msviper$regulons)[1:2] <- c("source", "target")

NR_R_Good_n_Bad_aracn_msviper_network <- Make_specific_network(Net = NR_R_Good_n_Bad_aracn_msviper$regulons, DEGs = RNAseq_diff_gene_expression_analysis$R_OR_NR_B$`NR.B-R.B`, msvip = NR_R_Good_n_Bad_aracn_msviper, title = "Good_n_Bad_aracn_network", collection = "NR_vs_R", P_P = F, PPI = F)

plot(NR_R_Good_n_Bad_aracn_msviper$mrs, cex = .7)






regulons_Baseline <- Baseline_aracne_network
colnames(regulons_Baseline)[1:2] <- c("tf", "target")

NR_R_Baseline_aracn_msviper <- run_msviper(RNAseq, regulons_Baseline, use_aracne = T, ref_R_B, ref_NR_B,  "R", "NR", minsize = 4, ges.filter=T)
colnames(NR_R_Baseline_aracn_msviper$regulons)[1:2] <- c("source", "target")

NR_R_Baseline_aracn_msviper_network <- Make_specific_network(Net = NR_R_Baseline_aracn_msviper$regulons, DEGs = RNAseq_diff_gene_expression_analysis$R_OR_NR_B$`NR.B-R.B`, msvip = NR_R_Baseline_aracn_msviper, title = "Baseline_aracn_network", collection = "NR_vs_R", P_P = F, PPI = F)

plot(NR_R_Baseline_aracn_msviper$mrs, cex = .7)

mrshadow <- shadow(NR_R_Baseline_aracn_msviper$mrs, regulators = 25, verbose = FALSE)
plot(mrshadow)
summary(mrshadow)


regulonaml_msviper <- regulonaml_SYMBOL
colnames(regulonaml_msviper) <- c("tf", "target", "mor", "likelihood")

NR_R_C_confidence_msviper <- run_msviper(RNAseq, regulons, use_aracne = T, ref_R_B, ref_NR_B,  "R", "NR", minsize = 4, ges.filter=T)





createNetworkFromDataFrames(edges = Bad_responder_aracne_network, title = "Bad_responder", collection = "Aracne network")
createNetworkFromDataFrames(edges = Baseline_aracne_network, title = "Baseline", collection = "Aracne network")
createNetworkFromDataFrames(edges = Good_n_Bad_aracne_network, title = "Good_n_Bad", collection = "Aracne network")

