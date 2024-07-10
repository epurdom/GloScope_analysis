sce <- readRDS("../../data/Processed_Datasets/perezLupus/perezLupus_default/perezLupus_default.Rds")
library(SingleCellExperiment)
#load("../../results/Perez2022/batchid_harmony.Rda")
cell_sample_ids <- sce$sample
pca_embedding <- reducedDim(sce)
library(umap)
pca_umap <-umap(pca_embedding[,1:10])

save(pca_umap, file = "../../results/BatchStudy/perez_umap_update.Rda")