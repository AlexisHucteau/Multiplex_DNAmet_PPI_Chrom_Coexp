rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(chaser)
  library(stringr)
  library(biomaRt)
  library(gridExtra)
  library(grid)
  library(ggplot2)
})

setwd("~/GitHub/Multiplex_DNAmet_PPI_Chrom_Coexp/")
source("~/Core_scripts/core_functions.R")
source("~/GitHub/Multiplex_DNAmet_PPI_Chrom_Coexp/Make_koichi_factor.R")

################### FUNCTIONS ####################

prepare_pchic <- function(cell_lines = NULL, minimum_interaction = 5, pchic = NULL){
  if (is.null(pchic)){
    load("DATA/pchic.RData")
  }
  if (is.null(cell_lines)){
    cell_lines = c("Mon", "Mac0", "Mac1", "Mac2", "Neu", "MK", "EP", "Ery", "FoeT", "nCD4", "tCD4", "aCD4", "naCD4", "nCD8", "tCD8", "nB", "tB")
  }
  pchic <- data.frame(pchic[rowSums(pchic[,cell_lines] >= minimum_interaction) >= 1,c(1,2,3,6,7,8)]) %>% na.omit(.)
  return(pchic)
}

################## DATAS ###########################

RNAseq_1 <- read.csv("~/GitHub/Koichi_gene_expression_analyses_git/Koichi_gene_expression_analyses/DATA/RNAseq_parsed.csv", row.names = 1, header = T, check.names = F)
RNAseq <- read.csv("~/GitHub/Koichi_gene_expression_analyses_git/Koichi_gene_expression_analyses/DATA/RNAseq.csv", row.names = 1, header = T, check.names = F)
# 
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# gene_ensembl <- biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "hgnc_symbol", values = rownames(RNAseq_1), mart = ensembl)
# 
# gene_ensembl_2 <- gene_ensembl[!duplicated(gene_ensembl$hgnc_symbol) & gene_ensembl$hgnc_symbol != "CCL3L3",]
# 
# RNAseq_1_ensembl <- RNAseq_1[gene_ensembl_2$hgnc_symbol,]
# rownames(RNAseq_1_ensembl) <- gene_ensembl_2$ensembl_gene_id

cell_lines <- c("Mon", "Mac1", "Mac0", "Mac2", "MK", "Ery", "EP")

load("DATA/pchic.RData")
pchic$baitChr <- paste0("chr", pchic$baitChr)
pchic$oeChr <- paste0("chr", pchic$oeChr)

all_cell_lines <- colnames(pchic)[12:28]

pchics <- lapply(all_cell_lines, function(cells){
  ps <- prepare_pchic(cells, pchic = pchic)
  chaser::make_chromnet(ps)
})
names(pchics) <- all_cell_lines
# 
conv=read.delim(paste('/mnt/SERVER-CRCT-STORAGE/CRCT21/Private/common/forGARDEN-NET/biomart/','BLUEPRINT_fragments_good.tsv', sep='/'), sep='\t')
# 
convpro=data.frame(paste('chr',conv$chr, ':', conv$start, '-', conv$end, sep=''), conv$ensembl)
rownames(convpro)=convpro[,1]
colnames(convpro)[2]='ens.id'
# 
# t <- merge(convpro, RNAseq_1_ensembl, by.x="ens.id", by.y = 0, all.x = F, all.y = F)
# rownames(t) <- t[,2]
# t<-t[,-c(1,2)]
# t <- data.frame(t, rowMeans(t), check.names = F)

t <- read.csv("~/t.csv", row.names=1, check.names = F)

setwd("./Results/")




chas_for_selected_pheno <- function(pchic = pchics, RNAseq = t, pheno, title, cells){
  RNAseq <- RNAseq[,pheno]
  
  pchics_pp <- lapply(names(pchics), function(cell_type){
    message(cell_type)
    baits <- export(pchics[[cell_type]], 'baits')
    pp <- chaser::subset_chromnet(pchics[[cell_type]], method = "nodes", nodes1 = baits)
    pp_exp <- chaser::load_features(pp,RNAseq,type='features_on_nodes',featnames = colnames(RNAseq), missingv=0)
    chas <- chaser::chas(pp_exp)
    chas_random <- tryCatch(chaser::randomize(pp_exp, nrandom = 10, dist.match = T), error=function(e) NULL)
    chas_random <- lapply(chas_random, chaser::chas)
    feat <- chaser::export(pp_exp)
    list("pp_exp" = pp_exp, "chas" = chas, "rand" = chas_random, "feat" = feat)
  })
  names(pchics_pp) <- cells
  chas_cell_lines <- sapply(names(pchics_pp), function(cell_type){
    pchics_pp[[cell_type]][["chas"]]
  })
  message("Pchic_pp DONE!")
  chas_cell_lines_df <- data.frame(cell_lines = rep(cells, each = nrow(chas_cell_lines)),
                                   values = as.vector(chas_cell_lines))
  
  p <- ggplot2::ggplot(chas_cell_lines_df, aes(x = cell_lines, y = values, fill = cell_lines))
  p <- p + ggplot2::geom_violin()
  p <- p + theme(axis.line = element_line(colour = "black"),
                 axis.text.x=element_text(size=16),
                 axis.text.y=element_text(size=16),
                 axis.title.x=element_text(size=16),
                 axis.title.y=element_text(size=16))
  p <- p + ggplot2::geom_jitter()
  p <- p + ggtitle(title)
  
  list("plot" = p,
       "pchics_pp" = pchics_pp)
}




Factor_R_OR_NR_B <- stringr::str_replace(Factor_R_OR_NR_B, pattern = "OR.REL", replacement = "REL")
Factor_R_OR_NR_B <- stringr::str_replace(Factor_R_OR_NR_B, pattern = "R.REL", replacement = "REL")
Factor_R_OR_NR_B <- as.factor(Factor_R_OR_NR_B)

Pheno_chas <- lapply(levels(Factor_R_OR_NR_B), function(pheno){
  ref <- Factor_R_OR_NR_B == pheno
  chas_for_selected_pheno(pheno = ref, title = pheno, cells = all_cell_lines)
})
names(Pheno_chas) <- levels(Factor_R_OR_NR_B)






png(file = paste0("./Chas_GE_all_Clines.png"), width = 1920, height = 1080)
grid.arrange(Pheno_chas[["NR.B"]]$plot, Pheno_chas[["R.B"]]$plot, Pheno_chas[["REL"]]$plot,
             ncol=3,
             top = textGrob('Chas based on response',
                            just = c('center'),
                            gp = gpar(fontsize = 32)))
dev.off()




Rands <- lapply(names(Pheno_chas), function(Comp){
  tmp <- lapply(names(Pheno_chas[[Comp]][["pchics_pp"]]), function(cell_line){
    tmp <- unlist(Pheno_chas[[Comp]][["pchics_pp"]][[cell_line]][["rand"]]) %>% unname()
    data.frame(cell_lines = rep(cell_line, length(tmp)), 
                              values = tmp)
  })
  dplyr::bind_rows(tmp)
})
names(Rands) <- levels(Factor_R_OR_NR_B)





prand <- lapply(names(Rands), function(comp){
  prand <- ggplot2::ggplot(Rands[[comp]], aes(x = cell_lines, y = values, fill = cell_lines))
  prand <- prand + ggplot2::geom_violin()
  prand <- prand + theme(axis.line = element_line(colour = "black"),
             axis.text.x=element_text(size=16),
             axis.text.y=element_text(size=16),
             axis.title.x=element_text(size=16),
             axis.title.y=element_text(size=16))
  prand <- prand + ggplot2::geom_jitter()
  prand <- prand + ggtitle(comp)
})
names(prand) <- levels(Factor_R_OR_NR_B)







png(file = paste0("./Chas_GE_all_Clines_rand.png"), width = 1920, height = 1080)
grid.arrange(prand[["NR.B"]], prand[["R.B"]], prand[["REL"]],
  ncol=3,
  top = textGrob('Chas rand based on response',
  just = c('center'),
  gp = gpar(fontsize = 32)))
dev.off()







Combined_plot <- lapply(names(Pheno_chas), function(comp){
  n_samples <- unname(table(Pheno_chas[[comp]][["plot"]][["data"]]$cell_lines)[1])
  p <- ggplot2::ggplot(Pheno_chas[[comp]][["plot"]][["data"]], aes(x = cell_lines, y = values, fill = cell_lines)) 
  p <- p + ggplot2::geom_violin()
  p <- p + ggplot2::geom_violin(Rands[[comp]], inherit.aes = F, mapping = aes(x = cell_lines, y = values, fill = cell_lines))
  p <- p + ggplot2::geom_jitter(Rands[[comp]], inherit.aes = F, mapping = aes(x = cell_lines, y = values, fill = cell_lines))
  p <- p + theme(axis.line = element_line(colour = "black"),
                 axis.text.x=element_text(size=16),
                 axis.text.y=element_text(size=16),
                 axis.title.x=element_text(size=16),
                 axis.title.y=element_text(size=16))
  p <- p + ggtitle(paste0(comp, " n=", n_samples))
  p
})
names(Combined_plot) <- levels(Factor_R_OR_NR_B)

png(file = paste0("./Combined_Chas_all_Clines_GE.png"), width = 1920, height = 1080)
grid.arrange(Combined_plot[["NR.B"]], Combined_plot[["R.B"]], Combined_plot[["REL"]],
             ncol=3,
             top = textGrob('Chas based on response',
                            just = c('center'),
                            gp = gpar(fontsize = 32)))
dev.off()









################### zscore

Chas_rand_mean <- lapply(Rands, function(comp){
  tmp <- split(comp, f = comp$cell_lines)
  tmp <- sapply(tmp, function(cell_line){
    mean(cell_line$values)
  }) %>% as.data.frame() %>% t()
  tmp
})

Chas_rand_sd <- lapply(Rands, function(comp){
  tmp <- split(comp, f = comp$cell_lines)
  tmp <- sapply(tmp, function(cell_line){
    sd(cell_line$values)
  }) %>% as.data.frame() %>% t()
  tmp
})

Chas_GE_Myelo_zscore <- lapply(names(Pheno_chas), function(comp){
  tmp <- Pheno_chas[[comp]][["plot"]][["data"]]
  tmp <- split(tmp, tmp$cell_lines)
  cell_lines <- colnames(Chas_rand_mean[[comp]]) %>% .[which(. %in% c("Mon", "Mac1", "Mac0", "Mac2", "MK", "Ery", "EP"))]
  n <- nrow(tmp[[cell_lines[1]]])
  tmp <- lapply(cell_lines, function(cell_line){
    message(cell_line)
    tmp[[cell_line]] %>% mutate(zscore = (values - Chas_rand_mean[[comp]][,cell_line])/Chas_rand_sd[[comp]][,cell_line])
  })
  names(tmp) <- cell_lines
  print(tmp)
  tmp <- sapply(tmp, function(cell_line){
    as.vector(cell_line$zscore)
  })
  res <- data.frame(cell_lines = rep(cell_lines, each = n),
                    zscore = as.vector(tmp))
})
names(Chas_GE_Myelo_zscore) <- levels(Factor_R_OR_NR_B)






Chas_GE_Myelo_zscore_plot <- lapply(names(Chas_GE_Myelo_zscore), function(comp){
  n_samples <- unname(table(Chas_GE_Myelo_zscore[[comp]]$cell_lines)[1])
  p <- ggplot2::ggplot(Chas_GE_Myelo_zscore[[comp]], aes(x = cell_lines, y = zscore, fill = cell_lines)) 
  p <- p + ggplot2::geom_violin()
  p <- p + ggplot2::geom_jitter()
  p <- p + theme(axis.line = element_line(colour = "black"),
                 axis.text.x=element_text(size=16),
                 axis.text.y=element_text(size=16),
                 axis.title.x=element_text(size=16),
                 axis.title.y=element_text(size=16))
  p <- p + ggtitle(paste0(comp, " n=", n_samples))
  p
})
names(Chas_GE_Myelo_zscore_plot) <- levels(Factor_R_OR_NR_B)





png(file = paste0("./Zscore_Chas_Myeloid_Clines_GE.png"), width = 1920, height = 1080)
grid.arrange(Chas_GE_Myelo_zscore_plot[["NR.B"]], Chas_GE_Myelo_zscore_plot[["R.B"]], Chas_GE_Myelo_zscore_plot[["REL"]],
             ncol=3,
             top = textGrob('Chas GE zscore based on response Myeloid cells',
                            just = c('center'),
                            gp = gpar(fontsize = 32)))
dev.off()



######### LSC17 score

LSC17_scores <- sapply(colnames(RNAseq_1), function(patient){
  DNMT3B <- RNAseq_1["DNMT3B", patient] * 0.0874
  ZBTB46 <- RNAseq_1["ZBTB46", patient] * -0.0347
  NYNRIN <- RNAseq_1["NYNRIN", patient] * 0.00865
  ARHGAP22 <- RNAseq_1["ARHGAP22", patient] * -0.0138
  LAPTM4B <- RNAseq_1["LAPTM4B", patient] * 0.00582
  MMRN1 <- RNAseq_1["MMRN1", patient] * 0.0258
  DPYSL3 <- RNAseq_1["DPYSL3", patient] * 0.0284
  KIAA0125 <- RNAseq_1["KIAA0125", patient] * 0.0196
  CDK6 <- RNAseq_1["CDK6", patient] * -0.0704
  CPXM1 <- RNAseq_1["CPXM1", patient] * -0.0258
  SOCS2 <- RNAseq_1["SOCS2", patient] * 0.0271
  EMP1 <- RNAseq_1["EMP1", patient] * 0.0146
  NGFRAP1 <- RNAseq_1["NGFRAP1", patient] * 0.0465
  CD34 <- RNAseq_1["CD34", patient] * 0.0338
  AKR1C3 <- RNAseq_1["AKR1C3", patient] * -0.0402
  GPR56 <- RNAseq_1["GPR56", patient] * 0.0501

  LSC17_score <- DNMT3B + ZBTB46 + NYNRIN +  ARHGAP22 + LAPTM4B + MMRN1 + DPYSL3 + KIAA0125 + CDK6 + CPXM1 + SOCS2 + EMP1 + NGFRAP1 + CD34 + AKR1C3 + GPR56
  LSC17_score
})


LSC17_scores_zscore <- (LSC17_scores - mean(LSC17_scores)) / sd(LSC17_scores)

LSC17_pheno <- ifelse(LSC17_scores_zscore >= 0, "LSC17_high", "LSC17_low")

All_dataset_chas <- chas_for_selected_pheno(pchic = pchics, RNAseq = t, pheno = rep(T, ncol(t)), title = "all dataset", cells = all_cell_lines)

Zscore_LSC17_chas <- lapply(names(All_dataset_chas$pchics_pp), function(cell_type){
  if(length(All_dataset_chas$pchics_pp[[cell_type]]$rand) == 0){
    return(NULL)
  }
  n <- seq(1:(length(All_dataset_chas$pchics_pp[[cell_type]]$rand[[1]]) - 1))
  rand_mean <- sapply(n, function(patient){
    sapply(All_dataset_chas$pchics_pp[[cell_type]]$rand, function(rand_chas){
      rand_chas[patient]
    }) %>% c() %>% mean()
  })
  rand_sd <- sapply(n, function(patient){
    sapply(All_dataset_chas$pchics_pp[[cell_type]]$rand, function(rand_chas){
      rand_chas[patient]
    }) %>% c() %>% sd()
  })
  
  Zscore = (All_dataset_chas$pchics_pp[[cell_type]]$chas[n] - rand_mean) / rand_sd

  data.frame(Patient_chas = All_dataset_chas$pchics_pp[[cell_type]]$chas[n],
             LSC17_score = LSC17_scores,
             Zscore = Zscore)
})

names(Zscore_LSC17_chas) <- names(pchics)


Zscore_LSC17_chas_plot <- lapply(names(Zscore_LSC17_chas), function(cell_type){
  p <- ggplot2::ggplot(Zscore_LSC17_chas[[cell_type]], aes(x = LSC17_score, y = Zscore))
  p <- p + ggplot2::geom_jitter()
  p <- p + theme(axis.line = element_line(colour = "black"),
                 axis.text.x=element_text(size=16),
                 axis.text.y=element_text(size=16),
                 axis.title.x=element_text(size=16),
                 axis.title.y=element_text(size=16))
  p <- p + ggtitle(cell_type)
  p
})
names(Zscore_LSC17_chas_plot) <- names(pchics)

LSC17_chas_plot <- lapply(names(Zscore_LSC17_chas), function(cell_type){
  p <- ggplot2::ggplot(Zscore_LSC17_chas[[cell_type]], aes(x = LSC17_score, y = Patient_chas))
  p <- p + ggplot2::geom_jitter()
  p <- p + theme(axis.line = element_line(colour = "black"),
                 axis.text.x=element_text(size=16),
                 axis.text.y=element_text(size=16),
                 axis.title.x=element_text(size=16),
                 axis.title.y=element_text(size=16))
  p <- p + ggtitle(cell_type)
  p
})
names(LSC17_chas_plot) <- names(pchics)

png(file = paste0("./Chas_all_cell_line_LSC17score_GE.png"), width = 1920, height = 1080)
grid.arrange(LSC17_chas_plot[["Mon"]], 
             LSC17_chas_plot[["Mac0"]], 
             LSC17_chas_plot[["Mac1"]],
             LSC17_chas_plot[["Mac2"]],
             LSC17_chas_plot[["EP"]],
             LSC17_chas_plot[["Ery"]],
             LSC17_chas_plot[["FoeT"]],
             LSC17_chas_plot[["nCD4"]],
             LSC17_chas_plot[["tCD4"]],
             LSC17_chas_plot[["naCD4"]],
             LSC17_chas_plot[["nCD8"]],
             LSC17_chas_plot[["nB"]],
             LSC17_chas_plot[["tB"]],
             LSC17_chas_plot[["Ery"]],
             ncol=4,
             top = textGrob('Chas GE based LSC17 score',
                            just = c('center'),
                            gp = gpar(fontsize = 32)))
dev.off()


png(file = paste0("./Zscore_Chas_all_cell_line_LSC17score_GE.png"), width = 1920, height = 1080)
grid.arrange(Zscore_LSC17_chas_plot[["Mon"]], 
             Zscore_LSC17_chas_plot[["Mac0"]], 
             Zscore_LSC17_chas_plot[["Mac1"]],
             Zscore_LSC17_chas_plot[["Mac2"]],
             Zscore_LSC17_chas_plot[["EP"]],
             Zscore_LSC17_chas_plot[["Ery"]],
             Zscore_LSC17_chas_plot[["FoeT"]],
             Zscore_LSC17_chas_plot[["nCD4"]],
             Zscore_LSC17_chas_plot[["tCD4"]],
             Zscore_LSC17_chas_plot[["naCD4"]],
             Zscore_LSC17_chas_plot[["nCD8"]],
             Zscore_LSC17_chas_plot[["nB"]],
             Zscore_LSC17_chas_plot[["tB"]],
             Zscore_LSC17_chas_plot[["Ery"]],
             ncol=3,
             top = textGrob('Chas GE zscore based LSC17 score',
                            just = c('center'),
                            gp = gpar(fontsize = 32)))
dev.off()









###################### Methylation

BMIQ <- readRDS("/media/alexis/DATA/Koichi_methylation_dat/BMIQ_norm_Koichi_samples.Rdata")
AnnoEpic <- read.csv("~/Illumina_Manifest/MethylationEPIC_v-1-0_B4.csv", skip = 7) %>% .[,c(2, 12, 13)]
Pheno <- read.csv("/media/alexis/DATA/Koichi_methylation_dat/samplesheet.csv") %>% unique()


AnnoBMIQ <- merge(BMIQ, AnnoEpic, by.x = 0, by.y = "Name")
AnnoBMIQ$end <- AnnoBMIQ$MAPINFO +1
AnnoBMIQ <- AnnoBMIQ[,c(107, 108, 109, 2:106)]
colnames(AnnoBMIQ)[1:3] <- c("chrom", "start", "end")
AnnoBMIQ$chrom <- paste0("chr", AnnoBMIQ$chrom)

write.table(AnnoBMIQ, "../DATA/Pchic/Meth_fragments.tsv", sep = "\t", row.names = F)




chas_met_for_selected_pheno <- function(pchic = pchics, BMIQ = AnnoBMIQ, pheno, title, cells){
  BMIQ <- BMIQ[,pheno]
  pchics_pp <- lapply(names(pchics), function(cell_type){
    message(cell_type)
    baits <- export(pchics[[cell_type]], 'baits')
    pp <- chaser::subset_chromnet(pchics[[cell_type]], method = "nodes", nodes1 = baits)
    pp_met <- chaser::load_features(pp,BMIQ,auxfun = 'mean', type='features_table',featnames = colnames(BMIQ), missingv=NA)
    pp_met<-subset_chromnet(pp_met, 'complete')
    chas <- chaser::chas(pp_met)
    chas_random <- tryCatch(chaser::randomize(pp_met, nrandom = 10, dist.match = T), error=function(e) NULL)
    chas_random <- lapply(chas_random, chaser::chas)
    feat <- chaser::export(pp_met)
    list("pp_met" = pp_met, "chas" = chas, "rand" = chas_random, "feat" = feat)
  })
  names(pchics_pp) <- cells
  chas_cell_lines <- sapply(names(pchics_pp), function(cell_type){
    pchics_pp[[cell_type]][["chas"]]
  })
  message("Pchic_pp DONE!")
  chas_cell_lines_df <- data.frame(cell_lines = rep(cells, each = nrow(chas_cell_lines)),
                                   values = as.vector(chas_cell_lines))
  
  p <- ggplot2::ggplot(chas_cell_lines_df, aes(x = cell_lines, y = values, fill = cell_lines))
  p <- p + ggplot2::geom_violin()
  p <- p + theme(axis.line = element_line(colour = "black"),
                 axis.text.x=element_text(size=16),
                 axis.text.y=element_text(size=16),
                 axis.title.x=element_text(size=16),
                 axis.title.y=element_text(size=16))
  p <- p + ggplot2::geom_jitter()
  p <- p + ggtitle(title)
  
  list("plot" = p,
       "pchics_pp" = pchics_pp)
}







Factor_R_OR_NR_B_met <- sapply(colnames(AnnoBMIQ)[4:108], function(phen){
  Pheno[which(Pheno$Sample == phen), "Pheno"][1]
}) 

Factor_R_OR_NR_B_met <- ifelse(stringr::str_detect(Factor_R_OR_NR_B_met, pattern = "Baseline.CR."), "R.B", 
              ifelse(stringr::str_detect(Factor_R_OR_NR_B_met, pattern = "Baseline.SD"), "NR.B", 
                     ifelse(stringr::str_detect(Factor_R_OR_NR_B_met, pattern = "Post_treatment"), "REL", "OR"))) %>% as.factor()







Pheno_met_chas <- lapply(levels(Factor_R_OR_NR_B_met), function(pheno){
  ref <- c(T, T, T)
  ref <- c(ref, Factor_R_OR_NR_B_met == pheno)
  ref <- c(ref[1:97], rep(F, 8))
  chas_met_for_selected_pheno(pheno = ref, title = pheno, cells = all_cell_lines)
})
names(Pheno_met_chas) <- levels(Factor_R_OR_NR_B_met)

png(file = paste0("./Chas_met_all_Clines.png"), width = 1920, height = 1080)
grid.arrange(Pheno_met_chas[["NR.B"]]$plot, Pheno_met_chas[["R.B"]]$plot, Pheno_met_chas[["REL"]]$plot,
             ncol=3,
             top = textGrob('Chas met based on response',
                            just = c('center'),
                            gp = gpar(fontsize = 32)))
dev.off()






Rands_met <- lapply(names(Pheno_met_chas), function(Comp){
  tmp <- lapply(names(Pheno_met_chas[[Comp]][["pchics_pp"]]), function(cell_line){
    tmp <- unlist(Pheno_met_chas[[Comp]][["pchics_pp"]][[cell_line]][["rand"]]) %>% unname()
    data.frame(cell_lines = rep(cell_line, length(tmp)), 
               values = tmp)
  })
  dplyr::bind_rows(tmp)
})

names(Rands_met) <- levels(Factor_R_OR_NR_B_met)






prand_met <- lapply(names(Rands_met), function(comp){
  prand <- ggplot2::ggplot(Rands_met[[comp]], aes(x = cell_lines, y = values, fill = cell_lines))
  prand <- prand + ggplot2::geom_violin()
  prand <- prand + theme(axis.line = element_line(colour = "black"),
                         axis.text.x=element_text(size=16),
                         axis.text.y=element_text(size=16),
                         axis.title.x=element_text(size=16),
                         axis.title.y=element_text(size=16))
  prand <- prand + ggplot2::geom_jitter()
  prand <- prand + ggtitle(comp)
})

names(prand_met) <- levels(Factor_R_OR_NR_B_met)







png(file = paste0("./Chas_met_all_Clines_rand.png"), width = 1920, height = 1080)
grid.arrange(prand_met[["NR.B"]], prand_met[["R.B"]], prand_met[["REL"]],
             ncol=3,
             top = textGrob('Chas rand met based on response',
                            just = c('center'),
                            gp = gpar(fontsize = 32)))
dev.off()






Chas_met_rand_sd <- lapply(Rands_met, function(comp){
  tmp <- split(comp, f = comp$cell_lines)
  tmp <- sapply(tmp, function(cell_line){
    sd(cell_line$values)
  }) %>% as.data.frame() %>% t()
  tmp
})

Chas_met_rand_mean <- lapply(Rands_met, function(comp){
  tmp <- split(comp, f = comp$cell_lines)
  tmp <- sapply(tmp, function(cell_line){
    mean(cell_line$values)
  }) %>% as.data.frame() %>% t()
  tmp
})

Chas_met_zscore <- lapply(names(Pheno_met_chas), function(comp){
  tmp <- Pheno_met_chas[[comp]][["plot"]][["data"]]
  tmp <- split(tmp, tmp$cell_lines)
  cell_lines <- colnames(Chas_met_rand_mean[[comp]]) %>% .[which(. %in% c("Mon", "Mac1", "Mac0", "Mac2", "MK", "Ery", "EP"))]
  n <- nrow(tmp[[cell_lines[1]]])
  tmp <- lapply(cell_lines, function(cell_line){
    message(cell_line)
    tmp[[cell_line]] %>% mutate(zscore = (values - Chas_met_rand_mean[[comp]][,cell_line])/Chas_met_rand_sd[[comp]][,cell_line])
  })
  names(tmp) <- cell_lines
  print(tmp)
  tmp <- sapply(tmp, function(cell_line){
    as.vector(cell_line$zscore)
  })
  res <- data.frame(cell_lines = rep(cell_lines, each = n),
                    zscore = as.vector(tmp))
})
names(Chas_met_zscore) <- levels(Factor_R_OR_NR_B_met)










Chas_met_zscore_plot <- lapply(names(Chas_met_zscore), function(comp){
  n_samples <- unname(table(Chas_met_zscore[[comp]]$cell_lines)[1])
  p <- ggplot2::ggplot(Chas_met_zscore[[comp]], aes(x = cell_lines, y = zscore, fill = cell_lines)) 
  p <- p + ggplot2::geom_violin()
  # p <- p + ggplot2::geom_jitter(Pheno_chas[[comp]][["plot"]][["data"]], inherit.aes = F, mapping =aes(x = cell_lines, y = values, fill = cell_lines))
  p <- p + theme(axis.line = element_line(colour = "black"),
                 axis.text.x=element_text(size=16),
                 axis.text.y=element_text(size=16),
                 axis.title.x=element_text(size=16),
                 axis.title.y=element_text(size=16))
  p <- p + ggtitle(paste0(comp, " n=", n_samples))
  p
})
names(Chas_met_zscore_plot) <- levels(Factor_R_OR_NR_B_met)





png(file = paste0("./Zscore_Chas_met_Myelo_Clines_PP.png"), width = 1920, height = 1080)
grid.arrange(Chas_met_zscore_plot[["NR.B"]], Chas_met_zscore_plot[["R.B"]], Chas_met_zscore_plot[["REL"]],
             ncol=3,
             top = textGrob('Chas zscore of methylation based on response',
                            just = c('center'),
                            gp = gpar(fontsize = 32)))
dev.off()


save.image("~/09.05.2022.RData")


######################### Cell lines

suppressPackageStartupMessages({
  library(hugene20sttranscriptcluster.db)
  library(oligo)
  library(data.table)
})

# Preparation of transcriptomes
Annotations <- data.frame(ACCNUM = sapply(contents(hugene20sttranscriptclusterACCNUM), paste, collapse = ", "), 
                          SYMBOL = sapply(contents(hugene20sttranscriptclusterSYMBOL), paste, collapse = ", "), 
                          DESC = sapply(contents(hugene20sttranscriptclusterGENENAME), paste, collapse = ", "),
                          ENSEMBL = sapply(contents(hugene20sttranscriptclusterENSEMBL), paste, collapse = ", "))
SDRF <- read.csv("~/GitHub/Epigenomic_integration/Epigenomic_integration/DATA_RNAseq_Lucille/Samplesheet.sdrf.csv")
celFiles <- list.celfiles("~/GitHub/Epigenomic_integration/Epigenomic_integration/DATA_RNAseq_Lucille", full.names = T)
Transcriptomes <- read.celfiles(celFiles) %>% rma() %>% oligo::exprs() %>% as.matrix()
Transcriptomes <- merge(Annotations, Transcriptomes, by.x = 0, by.y = 0, all.y = T)
Transcriptomes_cleaned <- Transcriptomes[which(Transcriptomes$ENSEMBL != "NA"),]
Transcriptomes_cleaned <- Transcriptomes_cleaned[which(Transcriptomes_cleaned$ENSEMBL %in% convpro$ens.id),]
Transcriptomes_cleaned <- Transcriptomes_cleaned[!duplicated(Transcriptomes_cleaned$ENSEMBL),]
rownames(Transcriptomes_cleaned) <- Transcriptomes_cleaned$ENSEMBL
AnnotCellLines <- Transcriptomes_cleaned[,c(1:5)]
Transcriptomes_cleaned <- Transcriptomes_cleaned[,c(6:29)]

t_cell_lines <- merge(convpro, Transcriptomes_cleaned, by.x="ens.id", by.y = 0, all.x = F, all.y = F)

rownames(t_cell_lines) <- t_cell_lines[,2]
t_cell_lines<-t_cell_lines[,-c(1,2)]
t_cell_lines <- data.frame(t_cell_lines, rowMeans(t_cell_lines), check.names = F)

t_cell_lines <- read.csv("~/t_cell_lines.csv", row.names=1, check.names = F)


Phenotype_cell_lines <- paste(SDRF$Characteristics.cell.line., SDRF$Characteristics.genotype., SDRF$Characteristics.treatment., sep = ".")
Phenotype_cell_lines <- as.factor(Phenotype_cell_lines)

target_pheno <- c("HL60.Mut.DMF", "HL60.Mut.AGI5198", "MOLM14.Mut.DMF", "MOLM14.Mut.AGI5198")

Pheno_chas_Cell_lines <- lapply(target_pheno, function(pheno){
  message("==================================================================")
  message("==================================================================")
  message(paste0("========== ", pheno, " =========="))
  ref <- Phenotype_cell_lines == pheno
  chas_for_selected_pheno(pheno = ref, title = pheno, cells = all_cell_lines, RNAseq = t_cell_lines)
})
names(Pheno_chas_Cell_lines) <- target_pheno



setwd("./Cell_lines/")


png(file = paste0("./Chas_GE_all_Clines_cell_lines.png"), width = 1920, height = 1080)
grid.arrange(Pheno_chas_Cell_lines[["HL60.Mut.DMF"]]$plot, Pheno_chas_Cell_lines[["HL60.Mut.AGI5198"]]$plot, Pheno_chas_Cell_lines[["MOLM14.Mut.DMF"]]$plot, Pheno_chas_Cell_lines[["MOLM14.Mut.AGI5198"]]$plot,
             ncol=3,
             top = textGrob('Chas based on cell line',
                            just = c('center'),
                            gp = gpar(fontsize = 32)))
dev.off()




Rands_cell_lines <- lapply(target_pheno, function(Comp){
  tmp <- lapply(names(Pheno_chas_Cell_lines[[Comp]][["pchics_pp"]]), function(cell_line){
    tmp <- unlist(Pheno_chas_Cell_lines[[Comp]][["pchics_pp"]][[cell_line]][["rand"]]) %>% unname()
    data.frame(cell_lines = rep(cell_line, length(tmp)), 
               values = tmp)
  })
  dplyr::bind_rows(tmp)
})
names(Rands_cell_lines) <- target_pheno





prand_cell_lines <- lapply(names(Rands_cell_lines), function(comp){
  prand <- ggplot2::ggplot(Rands_cell_lines[[comp]], aes(x = cell_lines, y = values, fill = cell_lines))
  prand <- prand + ggplot2::geom_violin()
  prand <- prand + theme(axis.line = element_line(colour = "black"),
                         axis.text.x=element_text(size=16),
                         axis.text.y=element_text(size=16),
                         axis.title.x=element_text(size=16),
                         axis.title.y=element_text(size=16))
  prand <- prand + ggplot2::geom_jitter()
  prand <- prand + ggtitle(comp)
})
names(prand_cell_lines) <- target_pheno







png(file = paste0("./Chas_GE_all_Clines_cell_lines_rand.png"), width = 1920, height = 1080)
grid.arrange(prand_cell_lines[["HL60.Mut.DMF"]], prand_cell_lines[["HL60.Mut.AGI5198"]], prand_cell_lines[["MOLM14.Mut.DMF"]], prand_cell_lines[["MOLM14.Mut.AGI5198"]],
             ncol=3,
             top = textGrob('Chas rand based on cell lines',
                            just = c('center'),
                            gp = gpar(fontsize = 32)))
dev.off()







Combined_plot_cell_lines <- lapply(target_pheno, function(comp){
  n_samples <- unname(table(Pheno_chas_Cell_lines[[comp]][["plot"]][["data"]]$cell_lines)[1])
  p <- ggplot2::ggplot(Pheno_chas_Cell_lines[[comp]][["plot"]][["data"]], aes(x = cell_lines, y = values, fill = cell_lines)) 
  p <- p + ggplot2::geom_violin()
  p <- p + ggplot2::geom_violin(Rands_cell_lines[[comp]], inherit.aes = F, mapping = aes(x = cell_lines, y = values, fill = cell_lines))
  p <- p + ggplot2::geom_jitter(Rands_cell_lines[[comp]], inherit.aes = F, mapping = aes(x = cell_lines, y = values, fill = cell_lines))
  p <- p + theme(axis.line = element_line(colour = "black"),
                 axis.text.x=element_text(size=16),
                 axis.text.y=element_text(size=16),
                 axis.title.x=element_text(size=16),
                 axis.title.y=element_text(size=16))
  p <- p + ggtitle(paste0(comp, " n=", n_samples))
  p
})
names(Combined_plot_cell_lines) <- target_pheno

png(file = paste0("./Combined_Chas_all_Clines__cell_lines_GE.png"), width = 1920, height = 1080)
grid.arrange(Combined_plot_cell_lines[["HL60.Mut.DMF"]], Combined_plot_cell_lines[["HL60.Mut.AGI5198"]], Combined_plot_cell_lines[["MOLM14.Mut.DMF"]], Combined_plot_cell_lines[["MOLM14.Mut.AGI5198"]],
             ncol=3,
             top = textGrob('Chas based on response',
                            just = c('center'),
                            gp = gpar(fontsize = 32)))
dev.off()









################### zscore

Chas_rand_mean_cell_lines <- lapply(Rands_cell_lines, function(comp){
  tmp <- split(comp, f = comp$cell_lines)
  tmp <- sapply(tmp, function(cell_line){
    mean(cell_line$values)
  }) %>% as.data.frame() %>% t()
  tmp
})

Chas_rand_sd_cell_lines <- lapply(Rands_cell_lines, function(comp){
  tmp <- split(comp, f = comp$cell_lines)
  tmp <- sapply(tmp, function(cell_line){
    sd(cell_line$values)
  }) %>% as.data.frame() %>% t()
  tmp
})

Chas_GE_Myelo_cell_lines_zscore <- lapply(target_pheno, function(comp){
  tmp <- Pheno_chas_Cell_lines[[comp]][["plot"]][["data"]]
  tmp <- split(tmp, tmp$cell_lines)
  cell_lines <- colnames(Chas_rand_mean_cell_lines[[comp]]) %>% .[which(. %in% c("Mon", "Mac1", "Mac0", "Mac2", "MK", "Ery", "EP"))]
  n <- nrow(tmp[[cell_lines[1]]])
  tmp <- lapply(cell_lines, function(cell_line){
    message(cell_line)
    tmp[[cell_line]] %>% mutate(zscore = (values - Chas_rand_mean_cell_lines[[comp]][,cell_line])/Chas_rand_sd_cell_lines[[comp]][,cell_line])
  })
  names(tmp) <- cell_lines
  print(tmp)
  tmp <- sapply(tmp, function(cell_line){
    as.vector(cell_line$zscore)
  })
  res <- data.frame(cell_lines = rep(cell_lines, each = n),
                    zscore = as.vector(tmp))
})
names(Chas_GE_Myelo_cell_lines_zscore) <- target_pheno






Chas_GE_Myelo_cell_lines_zscore_plot <- lapply(names(Chas_GE_Myelo_cell_lines_zscore), function(comp){
  n_samples <- unname(table(Chas_GE_Myelo_cell_lines_zscore[[comp]]$cell_lines)[1])
  p <- ggplot2::ggplot(Chas_GE_Myelo_cell_lines_zscore[[comp]], aes(x = cell_lines, y = zscore, colour = cell_lines, size = 5)) 
  p <- p + ggplot2::geom_jitter()
  p <- p + theme(axis.line = element_line(colour = "black"),
                 axis.text.x=element_text(size=16),
                 axis.text.y=element_text(size=16),
                 axis.title.x=element_text(size=16),
                 axis.title.y=element_text(size=16))
  p <- p + ggtitle(paste0(comp, " n=", n_samples))
  p
})
names(Chas_GE_Myelo_cell_lines_zscore_plot) <- target_pheno





png(file = paste0("./Zscore_Chas_Myeloid_Clines_cell_lines_GE.png"), width = 1920, height = 1080)
grid.arrange(Chas_GE_Myelo_cell_lines_zscore_plot[["HL60.Mut.DMF"]], Chas_GE_Myelo_cell_lines_zscore_plot[["HL60.Mut.AGI5198"]], Chas_GE_Myelo_cell_lines_zscore_plot[["MOLM14.Mut.DMF"]], Chas_GE_Myelo_cell_lines_zscore_plot[["MOLM14.Mut.AGI5198"]],
             ncol=3,
             top = textGrob('Chas GE zscore based on cell line (Myeloid cells)',
                            just = c('center'),
                            gp = gpar(fontsize = 32)))
dev.off()
