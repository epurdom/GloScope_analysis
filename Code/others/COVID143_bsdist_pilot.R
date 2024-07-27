library(vegan)
library(dplyr)

load("../../results/BatchStudy/Covid143_pilot.Rda")
load("../../results/COVID_143/meta.Rda")
plot_df$Status = as.character(plot_df$Status)
plot_df$sample_id = as.character(plot_df$sample_id)
plot_df$Site = as.character(plot_df$Site)
meta <- unique(plot_df[, c("sample_id", "Status", "Site")])

sub_id <- meta$sample_id[meta$Status %in% c("Covid", "Healthy")]
meta_sub <- meta[meta$sample_id %in% sub_id,]

pilot_mat <- pilot_mat[sub_id, sub_id]
gloscope_proportion_mat <- gloscope_proportion_mat[sub_id, sub_id]
pilot_mat_denovo <- pilot_mat_denovo[sub_id, sub_id]
gloscope_proportion_mat_denovo <- gloscope_proportion_mat_denovo[sub_id, sub_id]

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
  boot_id <- unlist(get_bootstrap_ids(meta_sub, "Site", "Status"))
  meta_boot <- meta_sub[boot_id,]
  pilot <-  get_stat(pilot_mat[boot_id, boot_id], meta_boot, batch = "Site", condition= "Status")

  
  clusprop <-  get_stat(gloscope_proportion_mat[boot_id, boot_id], meta_boot, batch = "Site", condition= "Status")
  
  pilot_denovo <-  get_stat(pilot_mat_denovo[boot_id, boot_id], meta_boot, batch = "Site", condition= "Status")

  
  clusprop_denovo <-  get_stat(gloscope_proportion_mat_denovo[boot_id, boot_id], meta_boot, batch = "Site", condition= "Status")

  stat_list = rbind(stat_list,c(pilot,clusprop, pilot_denovo,clusprop_denovo))
}

save(stat_list, file = "../../results/BatchStudy/COVID143_BSdist_pilot.Rda")
