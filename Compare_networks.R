library(dplyr)
setwd("GitHub/Multiplex_DNAmet_PPI_Chrom_Coexp/")

# Aracne-AP networks

Bad_responder_aracne_network <- read.csv("output_Bad_responder/network.txt", sep = "\t") 
Good_responder_aracne_network <- read.csv("output_Good_responder/network.txt", sep = "\t")
Baseline_aracne_network <- read.csv("output_Baseline/network.txt", sep = "\t")

colnames(Bad_responder_aracne_network)[c(1,2)] <- c("source", "target")
colnames(Good_responder_aracne_network)[c(1,2)] <- c("source", "target")
colnames(Baseline_aracne_network)[c(1,2)] <- c("source", "target")

# msviper aracn networks

NR_R_msviper <- readRDS("~/GitHub/Koichi_gene_expression_analyses_git/Koichi_gene_expression_analyses/Results/Tables/NR_R_msviper.rds")

library(viper)
library(aracne.networks)
library(GeneAnswers)
library(org.Hs.eg.db)

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
