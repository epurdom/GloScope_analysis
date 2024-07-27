devtools::load_all(here::here("..","popPackage"))
library(dplyr)
library(MASS)
library(mclust)
#load("../../results/BatchStudy/fabreLiver_plot_df.Rda")
seurat_object <- readRDS("../../data/Processed_Datasets/fabreLiver/fabreLiver_default/fabreLiver_default.Rds")
plot_df <- cbind(seurat_object@meta.data,
      seurat_object@reductions$pca@cell.embeddings,
      seurat_object@reductions$harmony@cell.embeddings,
      seurat_object@reductions$harPaper@cell.embeddings)
rm(seurat_object)
cell_sample_ids <- as.character(plot_df$sample)



pca_embedding <- plot_df[, paste0("PC_", 1:10)]
harmony_embedding <- plot_df[, paste0("harmony_", 1:10)]
harmony_paper_embedding <- plot_df[, paste0("harPaper_", 1:10)]


dist_mat_GMM_PCA = gloscope(pca_embedding, cell_sample_ids,
                                   dens = "GMM", dist_mat = "KL",
                                   r = 10000, 
                                   BPPARAM= BiocParallel::MulticoreParam(2,RNGseed = 1))

dist_mat_KNN_PCA = gloscope(pca_embedding, cell_sample_ids,
                                   dens = "KNN", dist_mat = "KL",
                                   r = 10000, 
                                   BPPARAM= BiocParallel::MulticoreParam(2,RNGseed = 1))


dist_mat_GMM_har_sample = gloscope(harmony_embedding, cell_sample_ids,
                                    dens = "GMM", dist_mat = "KL",
                                    r = 10000, 
                                    BPPARAM= BiocParallel::MulticoreParam(2,RNGseed = 1))


dist_mat_KNN_har_sample = gloscope(harmony_embedding, cell_sample_ids,
                                   dens = "KNN", dist_mat = "KL",
                                   r = 10000, 
                                   BPPARAM= BiocParallel::MulticoreParam(2,RNGseed = 1))


dist_mat_GMM_har_paper = gloscope(harmony_paper_embedding, cell_sample_ids,
                                   dens = "GMM", dist_mat = "KL",
                                   r = 10000, 
                                   BPPARAM= BiocParallel::MulticoreParam(2,RNGseed = 1))


dist_mat_KNN_har_paper = gloscope(harmony_paper_embedding, cell_sample_ids,
                                   dens = "KNN", dist_mat = "KL",
                                   r = 10000, 
                                   BPPARAM= BiocParallel::MulticoreParam(2,RNGseed = 1))



save(dist_mat_GMM_PCA, dist_mat_KNN_PCA,
     dist_mat_GMM_har_sample, dist_mat_KNN_har_sample,
     dist_mat_GMM_har_paper, dist_mat_KNN_har_paper,
     file = "../../results/BatchStudy/dist_mat_fabre_liver_update.Rda")
