library(vegan)
library(dplyr)

load("../../results/BatchStudy/Perez_pilot.Rda")
load("../../results/BatchStudy/Perez2022_subgroup_meta.Rda")
meta$sample <- paste0(substr(meta$batch_id,1,nchar(meta$batch_id)-3), ".",meta$Processing_Cohort)

#pilot_mat <- pilot_mat[meta$sample, meta$sample]
gloscope_proportion_mat <- gloscope_proportion_mat[meta$sample, meta$sample]
gloscope_proportion_mat_denovo <- gloscope_proportion_mat_denovo[meta$sample, meta$sample]


O2 <- function(PERMANOVA){
  sumsq <- PERMANOVA$SumOfSqs
  n <- PERMANOVA$Df
  o2 <- (sumsq[1]-n[1]*(sumsq[2]/n[2]))/(sumsq[3] + sumsq[2]/n[2])
  return(o2)
}

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


get_stat <- function(distmat, meta, batch, condition){
  ANS_batch <- anosim(distmat, meta[,batch])
  ANS_condition <- anosim(distmat, meta[,condition])
  ANO_batch <- adonis2(distmat ~meta[,batch])
  ANO_condition <- adonis2(distmat ~meta[,condition])
  O2_batch = O2(ANO_batch)
  O2_condition = O2(ANO_condition)
  
  return(c(O2_batch, O2_condition, ANS_batch$statistic, ANS_condition$statistic))
}
stat_list = c()
for(i in 1:100){
  print(i)
  set.seed(i)
  boot_id <- unlist(get_bootstrap_ids(meta, "Processing_Cohort", "disease_state"))
  meta_boot <- meta[boot_id,]
  pilot <-  get_stat(pilot_mat[boot_id, boot_id], meta_boot, batch = "Processing_Cohort", condition= "disease_state")

  clusprop <-  get_stat(gloscope_proportion_mat[boot_id, boot_id], meta_boot, batch = "Processing_Cohort", condition= "disease_state")
    
      pilot_denovo <-  get_stat(pilot_mat_denovo[boot_id, boot_id], meta_boot, batch = "Processing_Cohort", condition= "disease_state")

  clusprop_denovo <-  get_stat(gloscope_proportion_mat_denovo[boot_id, boot_id], meta_boot, batch = "Processing_Cohort", condition= "disease_state")
    
       set.seed(i)
   boot_id <- unlist(get_bootstrap_ids(meta,batch = "subgroup", condition = "disease_state"))
   meta_boot <- meta[boot_id,]
      pilot_new <-  get_stat(pilot_mat[boot_id, boot_id], meta_boot, batch = "subgroup", condition= "disease_state")

  clusprop_new <-  get_stat(gloscope_proportion_mat[boot_id, boot_id], meta_boot, batch = "subgroup", condition= "disease_state")
    
          pilot_denovo_new <-  get_stat(pilot_mat_denovo[boot_id, boot_id], meta_boot, batch = "subgroup", condition= "disease_state")

  clusprop_denovo_new <-  get_stat(gloscope_proportion_mat_denovo[boot_id, boot_id], meta_boot, batch = "subgroup", condition= "disease_state")
  

  stat_list = rbind(stat_list,c(pilot,clusprop,pilot_denovo,clusprop_denovo,
                                pilot_new,clusprop_new, pilot_denovo_new,clusprop_denovo_new))
}

save(stat_list, file = "../../results/BatchStudy/Perez2022_BSdist_pilot.Rda")
