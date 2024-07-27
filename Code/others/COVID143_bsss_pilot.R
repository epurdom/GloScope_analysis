library(vegan)
library(dplyr)
library(cluster)
load("../../results/BatchStudy/Covid143_pilot.Rda")
load("../../results/BatchStudy/COVID_143_meta.Rda")
plot_df$Status = as.character(plot_df$Status)
plot_df$sample_id = as.character(plot_df$sample)
plot_df$Site = as.character(plot_df$batch)
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
  boot_id <- unlist(get_bootstrap_ids(meta_sub, "Site", "Status"))
  meta_boot <- meta_sub[boot_id,]
  pilot <-  get_ss(meta_boot, pilot_mat[boot_id, boot_id], batch = "Site", group= "Status")

  
  clusprop <-  get_ss( meta_boot, gloscope_proportion_mat[boot_id, boot_id],batch = "Site", group= "Status")
  
  pilot_denovo <-  get_ss( meta_boot,pilot_mat_denovo[boot_id, boot_id], batch = "Site", group= "Status")

  
  clusprop_denovo <-  get_ss(meta_boot, gloscope_proportion_mat_denovo[boot_id, boot_id], batch = "Site", group= "Status")

  stat_list = rbind(stat_list,c(pilot,clusprop, pilot_denovo,clusprop_denovo))
}

save(stat_list, file = "../../results/BatchStudy/COVID143_bsss_pilot.Rda")
