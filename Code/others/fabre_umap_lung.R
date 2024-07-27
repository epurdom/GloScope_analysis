library(Seurat)
library(harmony)
library(dplyr)

seurat_object <- readRDS("../../data/Processed_Datasets/fabreLung/fabreLung_default/fabreLung_default.Rds")

seurat_object  = seurat_object %>%
  RunUMAP(reduction = "pca", dims = 1:10)
plot_df <- cbind(seurat_object@meta.data,
                seurat_object@reductions$pca@cell.embeddings,
                seurat_object@reductions$umap@cell.embeddings)

save(plot_df, file = "../../results/BatchStudy/fabre_umap_lung.Rda")
