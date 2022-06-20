library(tidyverse)
library(matrixTests)
library(data.table)
library(stringr)
library(ggrepel)
library(stringi)
library(rvest)


seurat_object<-TUH93_D14
metadata<-as.data.frame(read.csv("rxnMeta.csv"))
cluster<-2
fileTitle<-"effect_FLT3_neg"

reactions<-t(read.csv("TUH93_D14_scores.csv", row.names = 1, check.names = F))
reactions<-scale(reactions)
reactions<-cbind(high=NA,reactions)
reactions<-as.data.frame(reactions)
reactions$high <- ifelse(seurat_object@meta.data$seurat_clusters == cluster, 1, ifelse(seurat_object@meta.data$seurat_clusters == 1, 0, NA))

reactions<-as.data.frame(reactions)
names<-rownames(reactions)
reactions <- data.frame(apply(reactions, 2, function(x) as.numeric(as.character(x))))
rownames(reactions)<-names

test<-col_t_welch(subset(expFull, subset = day_200 == 1)[-1], subset(expFull, subset = day_200 == 0)[-1])
test.sig<-subset(test, subset = pvalue <= 0.05)

test.sig<-as.data.frame(test.sig)
test.sig<-t(test.sig)
colnames(test.sig)<-str_remove(colnames(test.sig), "_pos")
colnames(test.sig)<-str_remove(colnames(test.sig), "_neg")
#colnames(test.sig)<-str_remove(colnames(test.sig), "\\.e\\.")
for (i in 1:length(colnames(test.sig))){if(substring(colnames(test.sig)[i],1,1) == "X"){colnames(test.sig)[i]<-substring(colnames(test.sig)[i],2)}}
for (i in 1:length(colnames(test.sig))){if(substring(colnames(test.sig)[i],1,2) == "R_"){colnames(test.sig)[i]<-substring(colnames(test.sig)[i],3)}}
for (i in 1:length(colnames(test.sig))){if(stri_sub(colnames(test.sig)[i],-2) == "_e"){colnames(test.sig)[i]<-substring(colnames(test.sig)[i],1,stri_length(colnames(test.sig)[i])-2)}}

uni.test.sig <- test.sig[ , !duplicated(colnames(test.sig))]
uni.test.sig <- t(uni.test.sig)
names<-rownames(uni.test.sig)
uni.test.sig <- data.frame(apply(uni.test.sig, 2, function(x) as.numeric(as.character(x))))
rownames(uni.test.sig)<-names

name<-c()
subsystem<-c()
keggId<-c()

for (title in names){subsystem<-c(subsystem,metadata[metadata$id %like% title, ]$subsystem[1])}
for (title in names){name<-c(name, metadata[metadata$id %like% title, ]$name[1])}
for (title in names){keggId<-c(keggId, metadata[metadata$id %like% title, ]$keggId[1])}

fullNames<-paste0(name,keggId)
uni.test.sig$name<-fullNames
uni.test.sig$RECON3D_code<-rownames(uni.test.sig)

uni.test.sig$reaction_type<-subsystem
uni.test.sig$reaction_type[which(str_detect(rownames(uni.test.sig), "EX_"))] <- "EX/DEM"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$name, "endoplasmic reticular"))] <- "RT"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$name, "mitochondrial"))] <- "MT"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$name, "peroxisomal"))] <- "PT"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$name, "extracellular"))] <- "ET"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$name, "nuclear"))] <- "NT"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$name, "golgi apparatus"))] <- "GT"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$name, "lysosomal"))] <- "LT"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$name, " Aminosugars metabolism"))] <- "Aminosugars metabolism"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$name, " extracellular") & str_detect(uni.test.sig$name, "uniport"))] <- "ET -- UP"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$name, " extracellular") & str_detect(uni.test.sig$name, "MFS"))] <- "ET -- MF"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$name, " extracellular") & str_detect(uni.test.sig$name, "diffusion"))] <- "ET -- DIF"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$name, " extracellular") & str_detect(uni.test.sig$name, "sodium symport"))] <- "ET -- SS"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$name, " extracellular") & str_detect(uni.test.sig$name, "Organocation"))] <- "ET -- AAPO"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$name, " extracellular") & str_detect(uni.test.sig$name, "LAT"))] <- "ET -- LAT"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$reaction_type, "Citric acid cycle"))] <- "TCA"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$reaction_type, "Fatty acid oxidation"))] <- "FAO"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$reaction_type, "Fatty acid synthesis"))] <- "FAS"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$reaction_type, "Glycolysis"))] <- "GLYC"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$reaction_type, "Pentose phosphate pathway"))] <- "PPP"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$reaction_type, "Bile acid synthesis"))] <- "BAS"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$reaction_type, "Glutamate metabolism"))] <- "GLUT"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$reaction_type, "Nucleotide interconversion"))] <- "NI"

uni.test.sig$name<-str_replace(uni.test.sig$name," extracellular", "extracellular ")
uni.test.sig$name<-str_replace(uni.test.sig$name,"extracellular", "extracellular ")
uni.test.sig$name<-str_replace(uni.test.sig$name," endoplasmic reticular", "endoplasmic reticular ")
uni.test.sig$name<-str_replace(uni.test.sig$name,"endoplasmic reticular", "endoplasmic reticular ")
uni.test.sig$name<-str_replace(uni.test.sig$name," mitochondrial", "mitochondrial ")
uni.test.sig$name<-str_replace(uni.test.sig$name,"mitochondrial", "mitochondrial ")
uni.test.sig$name<-str_replace(uni.test.sig$name," nuclear", "nuclear ")
uni.test.sig$name<-str_replace(uni.test.sig$name," nuclear", "nuclear ")
uni.test.sig$name<-str_replace(uni.test.sig$name," golgi apparatus", "golgi apparatus ")
uni.test.sig$name<-str_replace(uni.test.sig$name,"golgi apparatus", " golgi apparatus")
uni.test.sig$name<-str_replace(uni.test.sig$name," peroxisomal", "peroxisomal ")
uni.test.sig$name<-str_replace(uni.test.sig$name,"peroxisomal ", "peroxisomal ")
uni.test.sig$name<-str_replace(uni.test.sig$name," lysosomal", "lysosomal ")
uni.test.sig$name<-str_replace(uni.test.sig$name,"lysosomal", "lysosomal ")

uni.test.sig<-uni.test.sig[order(uni.test.sig$statistic, decreasing = TRUE),]  

uni.test.sig$statistic<-as.numeric(uni.test.sig$statistic)
pos<-as.data.frame(head(uni.test.sig[which(uni.test.sig$statistic>0),],500))
neg<-as.data.frame(tail(uni.test.sig[which(uni.test.sig$statistic<0),],500))
#neg<- neg[seq(dim(neg)[1],1),]

web_scraper <-function(reaction_id){
  simple <- read_html(paste0("http://bigg.ucsd.edu/models/Recon3D/reactions/", reaction_id))
  scraped_data<-simple %>% html_nodes("p") %>% html_text()
  print(paste0("Original ID: ",reaction_id))
  print(paste0("Replacement: ",scraped_data[1]))
  print(paste0("Subsystem: ",scraped_data[6]))
  return(list(scraped_data[1],scraped_data[3],scraped_data[6]))
}

for (i in 1:length(rownames(pos))){if(str_detect(pos[i,18], "NANA")|str_detect(pos[i,18], "RE")| str_detect(pos[i,18], "HMR")& pos[i,19] != "age"){pos[i,18]<-ifelse(str_detect(pos[i,19], "EX_"),web_scraper(paste0(pos[i,19],"_e"))[1],web_scraper(pos[i,19])[1])
pos[i,20]<-ifelse(str_detect(pos[i,19], "EX_"),web_scraper(paste0(pos[i,19],"_e"))[2],web_scraper(pos[i,19])[2])}}
for (i in 1:length(rownames(neg))){if(str_detect(neg[i,18], "NANA")|str_detect(neg[i,18], "RE")| str_detect(neg[i,18], "HMR")& neg[i,19] != "age"){neg[i,18]<-ifelse(str_detect(neg[i,19], "EX_"),web_scraper(paste0(neg[i,19],"_e"))[1],web_scraper(neg[i,19])[1])
neg[i,20]<-ifelse(str_detect(neg[i,19], "EX_"),web_scraper(paste0(neg[i,19],"_e"))[2],web_scraper(neg[i,19])[2])}}

uni.test.sig<-rbind(pos, neg)

uni.test.sig$reaction_type[which(str_detect(rownames(uni.test.sig), "EX_"))] <- "EX/DEM"
uni.test.sig$reaction_type[which(str_detect(rownames(uni.test.sig), "SK_"))] <- "EX/DEM"
uni.test.sig$reaction_type[which(str_detect(rownames(uni.test.sig), "DM_"))] <- "EX/DEM"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$name, "endoplasmic reticular"))] <- "RT"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$name, "mitochondrial"))] <- "MT"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$name, "peroxisomal"))] <- "PT"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$name, "extracellular"))] <- "ET"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$reaction_type, "extracellular"))] <- "ET"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$reaction_type, "peroxisomal"))] <- "PT"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$name, "nuclear"))] <- "NT"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$name, "golgi apparatus"))] <- "GT"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$name, "lysosomal"))] <- "LT"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$reaction_type, "lysosomal"))] <- "LT"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$name, " Aminosugars metabolism"))] <- "Aminosugars metabolism"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$name, " extracellular") & str_detect(uni.test.sig$name, "uniport"))] <- "ET -- UP"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$name, " extracellular") & str_detect(uni.test.sig$name, "MFS"))] <- "ET -- MF"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$name, " extracellular") & str_detect(uni.test.sig$name, "diffusion"))] <- "ET -- DIF"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$name, " extracellular") & str_detect(uni.test.sig$name, "sodium symport"))] <- "ET -- SS"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$name, " extracellular") & str_detect(uni.test.sig$name, "Organocation"))] <- "ET -- AAPO"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$name, " extracellular") & str_detect(uni.test.sig$name, "LAT"))] <- "ET -- LAT"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$reaction_type, "Citric acid cycle"))] <- "TCA"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$reaction_type, "Fatty acid oxidation"))] <- "FAO"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$RECON3D_code, "FAOX"))] <- "FAO"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$reaction_type, "Fatty acid synthesis"))] <- "FAS"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$reaction_type, "Glycolysis"))] <- "GLYC"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$reaction_type, "Pentose phosphate pathway"))] <- "PPP"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$reaction_type, "Bile acid synthesis"))] <- "BAS"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$reaction_type, "Glutamate metabolism"))] <- "GLUT"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$reaction_type, "Nucleotide interconversion"))] <- "NI"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$RECON3D_code, "CHLST"))] <- "Cholesterol metabolism"
uni.test.sig$reaction_type[which(str_detect(uni.test.sig$RECON3D_code, "DESAT"))] <- "FAS"

write.csv(uni.test.sig, paste0(fileTitle,".csv"))
uni.test.sig<-read.csv(paste0(fileTitle,".csv"), header = T, row.names = 1)

pos<-as.data.frame(head(uni.test.sig[which(uni.test.sig$statistic>0),],100))
neg<-as.data.frame(tail(uni.test.sig[which(uni.test.sig$statistic<0),],100))

p <- ggplot(pos, aes(pos$statistic, -pos$pvalue, colour = reaction_type)) +  geom_point() +  
  geom_label_repel(aes(label = rownames(pos)),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50', max.overlaps = 20) +theme_classic()+ggtitle(paste0(fileTitle, " - Upregulated Reactions"))

q <- ggplot(neg, aes(neg$statistic, -neg$pvalue, colour = reaction_type)) +  geom_point() +  
  geom_label_repel(aes(label = rownames(neg)),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50', max.overlaps = 20) +theme_classic()+ggtitle(paste0(fileTitle, " - Downregulated Reactions"))


reactions_1<-rownames(subset(uni.test.sig, subset = statistic > 0))
reactions_1<-paste0("reaction:R_", reactions_1)
reactions_1<-str_c(reactions_1, collapse = "; ")

reactions_0<-rownames(subset(uni.test.sig, subset = statistic <= 0))
reactions_0<-paste0("reaction:R_", reactions_0)
reactions_0<-str_c(reactions_0, collapse = "; ")
