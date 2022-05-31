RNAseq <- read.csv("~/GitHub/Multiplex_DNAmet_PPI_Chrom_Coexp/DATA/GSE153348_IDH_RNA_Seq_matrix_submission.txt", sep = "\t", check.names = F)

gene_rownames <- rownames(RNAseq) %>% str_split(., "[|]") %>% lapply(., function(x){
  x[1]
}) %>%
  unlist(.) %>% 
  unique(.)


RNAseq$Gene <- rownames(RNAseq) %>% str_split(., "[|]") %>% lapply(., function(x){
  x[1]
}) %>%
  unlist(.)

RNAseq <- RNAseq %>%
  split(., .$Gene) %>%
  lapply(., function(x){
    l <- length(x[1,]) - 1
    cnames <- colnames(x)[c(1:l)]
    df <- x[,c(1:l)] %>%
      as.matrix(.) %>%
      colSums2(.) %>% 
      data.frame(.) %>%
      t(.) %>%
      data.frame(.)
    colnames(df) <- cnames
    df
  }) %>%
  rbindlist(.) %>% 
  data.frame(., check.names = F)

rownames(RNAseq) <- gene_rownames

RNAseq_CPM <- convertCounts(as.matrix(RNAseq),
                        unit = "CPM",
                        log = FALSE,
                        normalize = "none")

write.csv(RNAseq_CPM, "/mnt/SERVER-CRCT-STORAGE/CRCT18/UTILISATEURS/Alexis/DATA/RNAseq_mIDHi_cohort_CPM.csv")
