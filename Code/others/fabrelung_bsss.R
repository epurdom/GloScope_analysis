library(vegan)
library(dplyr)
library(stringr)
library(cluster)
load("../../results/BatchStudy/dist_mat_fabre_lung_update.Rda")
plot_df <- readRDS("../../data/Processed_Datasets/fabreLung/metadata/fabreLung_default_metadata.Rds")
meta <- unique(plot_df[, c("sample", "patient", "group", "batch")])



rownames(meta) <- as.character(meta$sample)
meta <- meta[rownames(dist_mat_GMM_PCA),]

meta$batch <- str_to_title(meta$batch)
meta$group <- str_to_title(meta$group)

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
  boot_id <- unlist(get_bootstrap_ids(meta, "batch", "group"))
  meta_boot <- meta[boot_id,]

    ss_pca <- get_ss(meta_boot, dist_mat_GMM_PCA[boot_id, boot_id], "batch", "group")
    ss_harsam <- get_ss(meta_boot, dist_mat_GMM_har_sample[boot_id, boot_id], "batch", "group")
    ss_harbat <- get_ss(meta_boot, dist_mat_GMM_har_paper[boot_id, boot_id], "batch", "group")
  
    ss_pca_knn <- get_ss(meta_boot, dist_mat_KNN_PCA[boot_id, boot_id], "batch", "group")
    ss_harsam_knn <- get_ss(meta_boot, dist_mat_KNN_har_sample[boot_id, boot_id], "batch", "group")
    ss_harbat_knn <- get_ss(meta_boot, dist_mat_KNN_har_paper[boot_id, boot_id], "batch", "group")

  stat_list = rbind(stat_list, c(ss_pca,ss_harsam, ss_harbat,
                                ss_pca_knn,ss_harsam_knn, ss_harbat_knn))
}

save(stat_list, file = "../../results/BatchStudy/fabreLung_BSss.Rda")
