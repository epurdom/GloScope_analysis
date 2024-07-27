library(vegan)
library(dplyr)
library(stringr)
library(cluster)

load("../../results/BatchStudy/Perez_pilot.Rda")
load("../../results/BatchStudy/Perez2022_subgroup_meta.Rda")
meta$sample <- paste0(substr(meta$batch_id,1,nchar(meta$batch_id)-3), ".",meta$Processing_Cohort)

#pilot_mat <- pilot_mat[meta$sample, meta$sample]
gloscope_proportion_mat <- gloscope_proportion_mat[meta$sample, meta$sample]
gloscope_proportion_mat_denovo <- gloscope_proportion_mat_denovo[meta$sample, meta$sample]



get_bootstrap_ids <- function(meta, batch, condition){
  bootstrap_sample_ids <- list()
  sample_ID <- 1:nrow(meta)
  
  for(i in unique(meta[,batch])){
    for(j in unique(meta[,condition])){
      # Generate bootstrap samples within the current class
      class_sample_ids <- sample_ID[meta[,batch] == i & meta[,condition] == j]
      
      bootstrap_sample_ids_cls <- sample(class_sample_ids, replace = TRUE)
      
      # Add the bootstrap sample IDs to the list
      bootstrap_sample_ids[[paste0(i,j)]] <- bootstrap_sample_ids_cls
      
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

stat_list = c()
for(i in 1:100){
  print(i)
  set.seed(i)
  boot_id <- unlist(get_bootstrap_ids(meta, "Processing_Cohort", "disease_state"))
  meta_boot <- meta[boot_id,]

  pilot <-  get_ss( meta_boot, pilot_mat[boot_id, boot_id],batch = "Processing_Cohort", group= "disease_state")
    

  clusprop <-  get_ss( meta_boot, gloscope_proportion_mat[boot_id, boot_id],batch = "Processing_Cohort", group= "disease_state")
    
      pilot_denovo <-  get_ss(meta_boot, pilot_mat_denovo[boot_id, boot_id], batch = "Processing_Cohort", group= "disease_state")

  clusprop_denovo <-  get_ss( meta_boot, gloscope_proportion_mat_denovo[boot_id, boot_id],batch = "Processing_Cohort", group= "disease_state")
    
       set.seed(i)
   boot_id <- unlist(get_bootstrap_ids(meta,batch = "subgroup", condition = "disease_state"))
   meta_boot <- meta[boot_id,]
      pilot_new <-  get_ss( meta_boot,pilot_mat[boot_id, boot_id], batch = "subgroup", group= "disease_state")

  clusprop_new <-  get_ss(meta_boot,gloscope_proportion_mat[boot_id, boot_id],  batch = "subgroup", group= "disease_state")
    
          pilot_denovo_new <-  get_ss( meta_boot, pilot_mat_denovo[boot_id, boot_id],batch = "subgroup", group= "disease_state")

  clusprop_denovo_new <-  get_ss(meta_boot, gloscope_proportion_mat_denovo[boot_id, boot_id], batch = "subgroup", group= "disease_state")
    
  stat_list = rbind(stat_list,c(pilot,clusprop,pilot_denovo,clusprop_denovo,
                                pilot_new,clusprop_new, pilot_denovo_new,clusprop_denovo_new))
}

save(stat_list, file = "../../results/BatchStudy/Perez2022_BSss_pilot.Rda")
