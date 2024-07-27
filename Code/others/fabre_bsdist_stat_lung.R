library(vegan)
library(dplyr)

load("../../results/BatchStudy/dist_mat_fabre_lung_update.Rda")
seurat_object <- readRDS("../../data/Processed_Datasets/fabreLung/fabreLung_default/fabreLung_default.Rds")
plot_df <- cbind(seurat_object@meta.data)
plot_df$sample <- as.character(plot_df$sample)
plot_df$batch <- as.character(plot_df$batch)
plot_df$group <- as.character(plot_df$group)

meta <- unique(plot_df[, c("sample", "group", "batch")])


# dist_mat_KNN_PCA[which(dist_mat_KNN_PCA<0)] = 0
# dist_mat_KNN_har_sample[which(dist_mat_KNN_har_sample<0)] = 0
# dist_mat_KNN_har_paper[which(dist_mat_KNN_har_paper<0)] = 0
# dist_mat_KNN_scvi[which(dist_mat_KNN_scvi<0)] = 0
# dist_mat_KNN_scvi_sample[which(dist_mat_KNN_scvi_sample<0)] = 0
# dist_mat_KNN_scvi_paper[which(dist_mat_KNN_scvi_paper<0)] = 0

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
  
  return(list(O2_batch, O2_condition, ANS_batch$statistic, ANS_condition$statistic))
}
stat_list = list()
for(i in 1:100){
  print(i)
  set.seed(i)
  boot_id <- unlist(get_bootstrap_ids(meta, "batch", "group"))
  meta_boot <- meta[boot_id,]
  GMM_PCA <- get_stat(dist_mat_GMM_PCA[boot_id, boot_id], meta_boot, batch = "batch", condition= "group")
  GMM_Har_sam <-  get_stat(dist_mat_GMM_har_sample[boot_id, boot_id], meta_boot, batch = "batch", condition= "group")
  GMM_Har_batch <-  get_stat(dist_mat_GMM_har_paper[boot_id, boot_id], meta_boot, batch = "batch", condition= "group")
  
  KNN_PCA <- get_stat(dist_mat_KNN_PCA[boot_id, boot_id], meta_boot, batch = "batch", condition= "group")
  KNN_Har_sam <-  get_stat(dist_mat_KNN_har_sample[boot_id, boot_id], meta_boot, batch = "batch", condition= "group")
  KNN_Har_batch <-  get_stat(dist_mat_KNN_har_paper[boot_id, boot_id], meta_boot, batch = "batch", condition= "group")

  
  

  stat_list[[i]] = list(GMM_PCA, GMM_Har_sam, GMM_Har_batch,
                        KNN_PCA, KNN_Har_sam, KNN_Har_batch)
}

save(stat_list, file = "../../results/BatchStudy/fabre_BSdist_stat_update.Rda")