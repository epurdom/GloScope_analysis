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
    stat_list = list(O2_batch, O2_condition, ANS_batch$statistic, ANS_condition$statistic)
  }
  
  return(stat_list)
}

 meta = readRDS("../../data/Processed_Datasets/stephensonCOVIDPBMC/metadata/stephensonCOVIDPBMC_default_metadata.Rds")
 meta$sample = as.character(meta$sample)
 meta$batch = as.character(meta$batch)
 
 meta <- unique(meta[, c("patient", "Status", "sample", "batch")])
 load("../../results/pseudobulk/pseudobulk_COVID143.Rda")
 plot_df$orig.ident = as.character(plot_df$orig.ident)
 plot_df$batch = meta$batch
 plot_df$group = meta$Status
 
 sub_id <- plot_df$orig.ident[plot_df$group %in% c("Covid", "Healthy")]
plot_df <- plot_df[plot_df$orig.ident %in% sub_id,]

 pseu_dist <- as.matrix(dist(plot_df[,c("PC_1","PC_2")]))
 
 load("../../results/BatchStudy/mofa_test_covid143.Rda")
 #meta$sample <- paste0(substr(meta$batch_id,1,nchar(meta$batch_id)-3), ".", meta$Processing_Cohort)
 all_factors <- MOFAcellulaR::get_tidy_factors(model = model,
                                               metadata = meta,
                                              factor = "all",
                                               sample_id_column = "sample")

 library(tidyr)
 
 
 data_wide <-  spread(all_factors, Factor, value)
 
 mofa_dist <- as.matrix(dist(data_wide[,5:9]))
 rownames(mofa_dist) <- colnames(mofa_dist) <- data_wide$sample

 sub_id <- which(data_wide$Status %in% c("Covid", "Healthy"))
 mofa_dist <- mofa_dist[sub_id, sub_id]
data_wide <- data_wide[sub_id,]
 
 
 mfa_meta <- as.data.frame(data_wide[,1:4])

stat_list = list()
for(i in 1:100){
  print(i)
  set.seed(i)
  boot_id <- unlist(get_bootstrap_ids(plot_df,batch = "batch", condition = "group"))
  meta_boot <- plot_df[boot_id,]
  pseudobulk <-  get_stat(pseu_dist[boot_id, boot_id], meta_boot, batch = "batch", condition= "group")
  
  set.seed(i)
  boot_id <- unlist(get_bootstrap_ids(mfa_meta,batch = "batch", condition = "Status"))
  meta_boot <- mfa_meta[boot_id,]
  MFA <-  get_stat(mofa_dist[boot_id, boot_id], meta_boot, batch = "batch", condition= "Status")


  
  stat_list[[i]] = list(pseudobulk,MFA)
}

save(stat_list, file = "../../results/BatchStudy/COVID143_bsdist_pseu_update.Rda")