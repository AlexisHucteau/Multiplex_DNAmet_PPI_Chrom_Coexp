Clinical_patient_data <- read.csv("~/GitHub/Koichi_gene_expression_analyses_git/Koichi_gene_expression_analyses/DATA/Clinical_patient_data.csv") %>%
  .[!duplicated(.),]

colRNA <- c("2292873_REL", "2410163_REL", "2297941_REL", "974025_REL", "2260553_BL", "2280935_BL", "2283265_BL", "2410163_BL", "2450523_BL", "2463247_BL", "2463247_REL", "2480485_BL", "2260671_REL", "2297625_BL", "2297625_REL", "1825001_BL", "1825001_REL", "2247707_BL", "2620771_BL", "2430935_BL", "2430935_REL", "2295595_BL", "2264225_BL", "2264225_REL", "2297941_BL", "2151503_BL", "2263055_BL", "757509_BL", "1794247_BL", "2010157_BL", "2150901_BL", "2152641_BL", "2163337_BL", "2210551_BL", "2219711_BL", "2292915_BL", "2296723_BL", "2361361_BL", "2361405_BL", "2366711_BL", "2370233_BL", "2382431_BL", "2419959_BL", "2498789_BL", "2502699_BL", "2512641_BL", "2583919_BL", "2616045_BL", "4128869_BL", "4235533_BL", "4259035_BL" )

Make_factor <- function(Samplesheet = Clinical_patient_data,
                        Samples_names,
                        met_data = F,
                        Mutations_to_ignore = 0,
                        Clinical_outcome_A = c("CR", "CRi"),
                        Clinical_name_A = "R",
                        Clinical_outcome_B = c("MLFS", "HI", "CRp", "PR"),
                        Clinical_name_B = "OR",
                        Clinical_outcome_C = c("SD", "PD"),
                        Clinical_name_C = "NR"){
  # Function made for Clinical_patient_data
  # Create a factor that can be used for Differential_analysis function
  # Samplesheet = Clinical_patient_data
  # Mutations_to_ignore: A vector of mutations that have to be taken into account (type 0 no mutations to ignore)
  # Clinical_outcome_A: A vector of best response corresponding to the phenotype A
  # Clinical_outcome_B: A vector of best response corresponding to the phenotype B
  # Clinical_outcome_C: A vector of best response corresponding to the phenotype C
  # Baseline_sample: A logical variable indicating whether Baseline samples are taken or not
  # Relapse_sample: A logical variable indicating whether Relapse samples are taken or not
  
  # Phenotype_A: The name of the first phenotype that have to be compared to
  # Phenotype _B: The name of the second phenotype that have to be compared to
  # Clinical_outcome_comparison: A logical variable indicating whether clinical outcome are taken into account
  # Baseline:
  # Relapse: A logical variable indicating whether Relapse samples are taken or not
  if(typeof(Mutations_to_ignore) != "double"){
    Mutations_samples <- Samplesheet[which(duplicated(str_split(Samplesheet$mutations, pattern=","), Mutations_to_ignore)),] %>%
      c(.$Baseline_RNAseq_data, .$Relapse_RNAseq_data) %>%
      na.omit()
    Mutations_factor <- factor(ifelse(Samples_names %in% Mutations_samples, "Mut", "WT"))
  }else{
    Mutations_factor <- factor(rep("", length(Samples_names)))
  }
  
  if(met_data){
    Clinical_outcome_A <- Samplesheet[which(Samplesheet$Best_response %in% Clinical_outcome_A),] %>%
      c(.$Baseline_Sample, .$Post_treatment_sample) %>%
      na.omit()
    Clinical_outcome_B <- Samplesheet[which(Samplesheet$Best_response %in% Clinical_outcome_B),] %>%
      c(.$Baseline_Sample, .$Post_treatment_sample) %>%
      na.omit()
    Clinical_outcome_C <- Samplesheet[which(Samplesheet$Best_response %in% Clinical_outcome_C),] %>%
      c(.$Baseline_Sample, .$Post_treatment_sample) %>%
      na.omit()
  }else{
    Clinical_outcome_A <- Samplesheet[which(Samplesheet$Best_response %in% Clinical_outcome_A),] %>%
      c(.$Baseline_RNAseq_data, .$Relapse_RNAseq_data) %>%
      na.omit()
    Clinical_outcome_B <- Samplesheet[which(Samplesheet$Best_response %in% Clinical_outcome_B),] %>%
      c(.$Baseline_RNAseq_data, .$Relapse_RNAseq_data) %>%
      na.omit()
    Clinical_outcome_C <- Samplesheet[which(Samplesheet$Best_response %in% Clinical_outcome_C),] %>%
      c(.$Baseline_RNAseq_data, .$Relapse_RNAseq_data) %>%
      na.omit()
  }
  Clinical_outcome <- factor(ifelse(Samples_names %in% Clinical_outcome_A, Clinical_name_A,
                                    ifelse(Samples_names %in% Clinical_outcome_B, Clinical_name_B,
                                           ifelse(Samples_names %in% Clinical_outcome_C, Clinical_name_C, ""))))
  if(met_data){
    Sample_timing <- factor(ifelse(Samples_names %in% Samplesheet$Baseline_Sample, "B", "Post"))
  }else{
    Sample_timing <- factor(ifelse(Samples_names %in% Samplesheet$Baseline_RNAseq_data, "B", "REL"))
  }
  if(typeof(Mutations_to_ignore) != "double"){
    Final_factor <- paste(Mutations_factor, Clinical_outcome, Sample_timing, sep = ".") %>% as.factor()
  }else{
    Final_factor <- paste(Clinical_outcome, Sample_timing, sep = ".") %>% as.factor()
  }
  
  return(Final_factor)
}


Factor_R_OR_NR_B <- Make_factor(Clinical_patient_data,
                                colRNA)

rm(colRNA)