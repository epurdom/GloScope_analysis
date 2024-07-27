library(stringr)
library(Seurat)
library(SingleCellExperiment)
library(muscat)

ddir = "../../data/Processed_Datasets/stephensonCOVIDPBMC/stephensonCOVIDPBMC_default/"

seurat_object <- readRDS(paste0(ddir, "stephensonCOVIDPBMC_default.Rds"))
seurat_object$cluster_id = as.character(seurat_object$cell_type)
seurat_object$sample_id = as.character(seurat_object$sample)

#pb_cluster <- aggregateData(sce, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))
#pb_samp <- aggregateData(sce, assay = "counts", fun = "sum", by = c("sample_id"))


#save(pb_cluster,pb_samp, file = "../../results/BatchStudy/Perez2022_pseudobulk.Rda")
pb_cluster <- list()
for(i in unique(seurat_object$cluster_id)){
  print(i)
  counts <- list()
  for(j in unique(seurat_object$sample_id)){
    ids <- seurat_object$cluster_id == i & seurat_object$sample_id == j
    if(sum(ids) == 0){
      counts[[j]] = 0
    }else if(sum(ids) == 1){
      tmp <- seurat_object@assays$raw@counts[,ids]
      counts[[j]] = tmp
    }else{
      tmp <- rowSums(seurat_object@assays$raw@counts[,ids])
      counts[[j]] = tmp
    }
  }
  counts <- do.call(cbind, counts)
  pb_cluster[[i]] = counts
}

pb_samp <- list()
for(i in unique(seurat_object$sample_id)){
    ids <-  seurat_object$sample_id == i
    tmp <- rowSums(seurat_object@assays$raw@counts[,ids])
    pb_samp[[i]] = tmp
}
pb_samp <- do.call(cbind, pb_samp)
save(pb_cluster,pb_samp, file = "../../results/BatchStudy/COVID143_pseudobulk_manual.Rda")

