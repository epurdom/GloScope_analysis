sil_score = c()
for(i in 1:100){
   load(paste0("../../results/simulation/tmp/sil/", i, ".Rda"))
#   print(i)
   plot_df = plot_df[plot_df$group_id == "A",]
   myclus = plot_df$pca
   load(paste0("../../../tmp/clusprop_V3/Setting1/25_50/",i, ".Rda"))
   plot_df = as.matrix(plot_df[plot_df$group_id == "A", paste0("PC_",1:10)])
   test = dist(plot_df)
   mysil = silhouette(as.numeric(myclus), test)
   sil_score <- c(sil_score, mean(mysil[,3]))
 }

save(sil_score, file = "../../results/simulation/cluster_sil/sil_pcaclus.Rda")


sil_score = c()
for(i in 1:100){
  load(paste0("../../results/simulation/tmp/sil/", i, ".Rda"))
  #   print(i)
  plot_df = plot_df[plot_df$group_id == "A",]
  myclus = plot_df$cluster_id
  load(paste0("../../../tmp/clusprop_V3/Setting1/25_50/",i, ".Rda"))
  plot_df = as.matrix(plot_df[plot_df$group_id == "A", paste0("PC_",1:10)])
  test = dist(plot_df)
  mysil = silhouette(as.numeric(myclus), test)
  sil_score <- c(sil_score, mean(mysil[,3]))
}

save(sil_score, file = "../../results/simulation/cluster_sil/sil_pcatrue.Rda")

sil_score = c()
for(i in 1:100){
  load(paste0("../../results/simulation/tmp/sil/", i, ".Rda"))
  #   print(i)
  plot_df = plot_df[plot_df$group_id == "A",]
  myclus = plot_df$sample_id
  load(paste0("../../../tmp/clusprop_V3/Setting1/25_50/",i, ".Rda"))
  plot_df = as.matrix(plot_df[plot_df$group_id == "A", paste0("PC_",1:10)])
  test = dist(plot_df)
  mysil = silhouette(as.numeric(myclus), test)
  sil_score <- c(sil_score, mean(mysil[,3]))
}

save(sil_score, file = "../../results/simulation/cluster_sil/sil_pcasam.Rda")

sil_score = c()
for(i in 1:100){
  load(paste0("../../results/simulation/tmp/sil/", i, ".Rda"))
  #   print(i)
  plot_df = plot_df[plot_df$group_id == "A",]
  myclus = plot_df$cluster_id
  load(paste0("../../../tmp/clusprop_V3/Setting1/25_50/",i, ".Rda"))
  plot_df = as.matrix(plot_df[plot_df$group_id == "A", paste0("scvi_",1:10)])
  test = dist(plot_df)
  mysil = silhouette(as.numeric(myclus), test)
  sil_score <- c(sil_score, mean(mysil[,3]))
}

save(sil_score, file = "../../results/simulation/cluster_sil/sil_scvitrue.Rda")

sil_score = c()
for(i in 1:100){
  load(paste0("../../results/simulation/tmp/sil/", i, ".Rda"))
  #   print(i)
  plot_df = plot_df[plot_df$group_id == "A",]
  myclus = plot_df$sample_id
  load(paste0("../../../tmp/clusprop_V3/Setting1/25_50/",i, ".Rda"))
  plot_df = as.matrix(plot_df[plot_df$group_id == "A", paste0("scvi_",1:10)])
  test = dist(plot_df)
  mysil = silhouette(as.numeric(myclus), test)
  sil_score <- c(sil_score, mean(mysil[,3]))
}

save(sil_score, file = "../../results/simulation/cluster_sil/sil_scvisam.Rda")

sil_score = c()
for(i in 1:100){
  load(paste0("../../results/simulation/tmp/sil/", i, ".Rda"))
  #   print(i)
  plot_df = plot_df[plot_df$group_id == "A",]
  myclus = plot_df$seurat_cluster
  load(paste0("../../../tmp/clusprop_V3/Setting1/25_50/",i, ".Rda"))
  plot_df = as.matrix(plot_df[plot_df$group_id == "A", paste0("scvi_",1:10)])
  test = dist(plot_df)
  mysil = silhouette(as.numeric(myclus), test)
  sil_score <- c(sil_score, mean(mysil[,3]))
}

save(sil_score, file = "../../results/simulation/cluster_sil/sil_scviclus.Rda")