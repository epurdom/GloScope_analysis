library(dplyr)
library(vegan)
library(stringr)
library(cluster)
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



get_ss <- function(meta, distmat, batch, group){
    ss_batch  <- silhouette(as.numeric(factor(meta[,batch])), distmat)
    ss_group  <- silhouette(as.numeric(factor(meta[,group])), distmat)
    
    return(c(batch = mean(ss_batch[,3]),
                group = mean(ss_group[,3])))
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

stat_list = c()
 for(i in 1:100){
   print(i)
   set.seed(i)
   boot_id <- unlist(get_bootstrap_ids(meta,batch = "Processing_Cohort", condition = "disease_state"))
   meta_boot <- meta[boot_id,]
   pseudobulk <-  get_ss(meta_boot, pseu_dist[boot_id, boot_id], batch = "Processing_Cohort", group= "disease_state")
   
   set.seed(i)
   boot_id <- unlist(get_bootstrap_ids(mfa_meta,batch = "Processing_Cohort", condition = "disease_state"))
   meta_boot <- mfa_meta[boot_id,]
   MFA <-  get_ss(meta_boot,mofa_dist[boot_id, boot_id],  batch = "Processing_Cohort", group= "disease_state")
 
   set.seed(i)
   boot_id <- unlist(get_bootstrap_ids(meta,batch = "subgroup", condition = "disease_state"))
   meta_boot <- meta[boot_id,]
   pseudobulk_new <-  get_ss(meta_boot, pseu_dist[boot_id, boot_id], batch = "subgroup", group= "disease_state")
   
   set.seed(i)
   boot_id <- unlist(get_bootstrap_ids(mfa_meta,batch = "subgroup", condition = "disease_state"))
   meta_boot <- mfa_meta[boot_id,]
   MFA_new <-  get_ss( meta_boot, mofa_dist[boot_id, boot_id],batch = "subgroup", group= "disease_state")
   
   
   stat_list = rbind(stat_list, c(pseudobulk,MFA,
                         pseudobulk_new,MFA_new))
 }

save(stat_list, file = "../../results/BatchStudy/Perez2022_bsss_pseu.Rda")