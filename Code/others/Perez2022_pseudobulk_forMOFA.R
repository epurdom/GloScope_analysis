library(stringr)
library(Seurat)
library(SingleCellExperiment)
library(muscat)

load("../../results/Perez2022/sce_object.Rda")
sce$sample_id = paste0(sce$sample_uuid,"_", sce$Processing_Cohort)
sce$cluster_id = as.character(sce$author_cell_type)

#pb_cluster <- aggregateData(sce, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))
#pb_samp <- aggregateData(sce, assay = "counts", fun = "sum", by = c("sample_id"))


#save(pb_cluster,pb_samp, file = "../../results/BatchStudy/Perez2022_pseudobulk.Rda")
pb_cluster <- list()
for(i in unique(sce$cluster_id)){
  print(i)
  counts <- list()
  for(j in unique(sce$sample_id)){
    ids <- sce$cluster_id == i & sce$sample_id == j
    if(sum(ids) == 0){
      counts[[j]] = 0
    }else if(sum(ids) == 1){
      tmp <- sce@assays@data$counts[,ids]
      counts[[j]] = tmp
    }else{
      tmp <- rowSums(sce@assays@data$counts[,ids])
      counts[[j]] = tmp
    }
  }
  counts <- do.call(cbind, counts)
  pb_cluster[[i]] = counts
}

pb_samp <- list()
for(i in unique(sce$sample_id)){
    ids <-  sce$sample_id == i
    tmp <- rowSums(sce@assays@data$counts[,ids])
    pb_samp[[i]] = tmp
}
pb_samp <- do.call(cbind, pb_samp)
save(pb_cluster,pb_samp, file = "../../results/BatchStudy/Perez2022_pseudobulk_manual.Rda")

