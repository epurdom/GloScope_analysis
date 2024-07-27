library(dplyr)
library(vegan)
library(stringr)
O2 <- function(PERMANOVA){
  sumsq <- PERMANOVA$SumOfSqs
  n <- PERMANOVA$Df
  o2 <- (sumsq[1]-n[1]*(sumsq[2]/n[2]))/(sumsq[3] + sumsq[2]/n[2])
  return(o2)
}

get_bootstrap_ids <- function(meta, batch = NULL, condition){
  bootstrap_sample_ids <- list()
  sample_ID <- 1:nrow(meta)
  
  if(is.null(batch)){
    for(j in unique(meta[,condition])){
      # Generate bootstrap samples within the current class
      class_sample_ids <- sample_ID[meta[,condition] == j]
      
      bootstrap_sample_ids_cls <- sample(class_sample_ids, replace = TRUE)
      
      # Add the bootstrap sample IDs to the list
      bootstrap_sample_ids[[j]] <- bootstrap_sample_ids_cls
    }
  }else{
    for(i in unique(meta[,batch])){
      for(j in unique(meta[,condition])){
        # Generate bootstrap samples within the current class
        class_sample_ids <- sample_ID[meta[,batch] == i & meta[,condition] == j]
        
        bootstrap_sample_ids_cls <- sample(class_sample_ids, replace = TRUE)
        
        # Add the bootstrap sample IDs to the list
        bootstrap_sample_ids[[paste0(i,j)]] <- bootstrap_sample_ids_cls
        
      }
    }
  }
  
  return(bootstrap_sample_ids)
}


get_stat <- function(distmat, meta, batch = NULL, condition){
  if(is.null(batch)){
    ANS_condition <- anosim(distmat, meta[,condition])
    ANO_condition <- adonis2(distmat ~meta[,condition])
    O2_condition = O2(ANO_condition)
    stat_list = list(O2_condition, ANS_condition$statistic)
    
  }else{
    ANS_batch <- anosim(distmat, meta[,batch])
    ANS_condition <- anosim(distmat, meta[,condition])
    ANO_batch <- adonis2(distmat ~meta[,batch])
    ANO_condition <- adonis2(distmat ~meta[,condition])
    O2_batch = O2(ANO_batch)
    O2_condition = O2(ANO_condition)
    stat_list = c(O2_batch, O2_condition, ANS_batch$statistic, ANS_condition$statistic)
  }
  
  return(stat_list)
}



meta = readRDS("../../data/Processed_Datasets/perezLupus/metadata/perezLupus_default_metadata.Rds")
meta <- unique(meta[, c("patient", "group", "sample", "batch")])
levels(meta$group)[which(levels(meta$group) == "na")] = "normal"

load("../../results/pseudobulk/pseudobulk_perez.Rda")

rownames(plot_df) = str_replace(rownames(plot_df), "_", ".")
meta <- meta[match(rownames(plot_df), meta$sample),]
plot_df = cbind(plot_df, meta)
plot_df$group = str_to_title(plot_df$group)

load("../../results/BatchStudy/Perez2022_subgroup_meta.Rda")
meta$sample <- paste0(substr(meta$batch_id,1,nchar(meta$batch_id)-3), ".",meta$Processing_Cohort)
rownames(meta) = meta$sample
meta <- meta[as.character(plot_df$sample),]
plot_df$subgroup <- meta$subgroup

pseu_dist <- as.matrix(dist(plot_df[,c("PC_1","PC_2")]))
 load("../../results/BatchStudy/mofa_test_perez.Rda")
 #meta$sample <- paste0(substr(meta$batch_id,1,nchar(meta$batch_id)-3), ".", meta$Processing_Cohort)
 all_factors <- MOFAcellulaR::get_tidy_factors(model = model,
                                               metadata = meta,
                                               factor = "all",
                                               sample_id_column = "sample")


 library(tidyr)
 data_wide <-  spread(all_factors, Factor, value)
 mofa_dist <- as.matrix(dist(data_wide[,7:11]))
 rownames(mofa_dist) <- colnames(mofa_dist) <- data_wide$sample
 mfa_meta <- as.data.frame(data_wide[,1:6])

 for(i in 1:100){
   print(i)
   set.seed(i)
   boot_id <- unlist(get_bootstrap_ids(meta,batch = "Processing_Cohort", condition = "disease_state"))
   meta_boot <- meta[boot_id,]
   pseudobulk <-  get_stat(pseu_dist[boot_id, boot_id], meta_boot, batch = "Processing_Cohort", condition= "disease_state")
   
   set.seed(i)
   boot_id <- unlist(get_bootstrap_ids(mfa_meta,batch = "Processing_Cohort", condition = "disease_state"))
   meta_boot <- mfa_meta[boot_id,]
   MFA <-  get_stat(mofa_dist[boot_id, boot_id], meta_boot, batch = "Processing_Cohort", condition= "disease_state")
 
   set.seed(i)
   boot_id <- unlist(get_bootstrap_ids(meta,batch = "subgroup", condition = "disease_state"))
   meta_boot <- meta[boot_id,]
   pseudobulk_new <-  get_stat(pseu_dist[boot_id, boot_id], meta_boot, batch = "subgroup", condition= "disease_state")
   
   set.seed(i)
   boot_id <- unlist(get_bootstrap_ids(mfa_meta,batch = "subgroup", condition = "disease_state"))
   meta_boot <- mfa_meta[boot_id,]
   MFA_new <-  get_stat(mofa_dist[boot_id, boot_id], meta_boot, batch = "subgroup", condition= "disease_state")
   
   
   stat_list[[i]] = list(pseudobulk,MFA,
                         pseudobulk_new,MFA_new)
 }

save(stat_list, file = "../../results/BatchStudy/Perez2022_bsdist_pseu.Rda")