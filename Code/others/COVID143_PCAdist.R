ddir = "../../results/COVID_143/"


library(Seurat)
library(dplyr)
library(MASS)
library(mclust)
library(BiocParallel)
library(GloScope)
#source("tmp_function.R")
seurat_object <- readRDS("../../data/Processed_Datasets/stephensonCOVIDPBMC/stephensonCOVIDPBMC_default/stephensonCOVIDPBMC_default.Rds")
plot_df <- cbind(seurat_object@meta.data,
                 seurat_object@reductions$pca@cell.embeddings,
                 seurat_object@reductions$pca_harmony@cell.embeddings)
save(plot_df, file = "../../results/BatchStudy/COVID_143_meta.Rda")
cell_sample_ids <- plot_df$sample
pca_embedding <- plot_df[, paste0("PC_", 1:10)]
harmony_embedding <- plot_df[,paste0("pca_harmony_",1:10)]

dist_mat_GMM_PCA = gloscope(pca_embedding, cell_sample_ids,
                            dens = "GMM", dist_mat = "KL",
                            r = 10000, 
                            BPPARAM= BiocParallel::SerialParam(RNGseed=2))


dist_mat_KNN_PCA = gloscope(pca_embedding, cell_sample_ids,
                            dens = "KNN", dist_mat = "KL",
                            r = 10000, k=50,
                            BPPARAM= BiocParallel::SerialParam(RNGseed=2))



dist_mat_GMM_har_sample = gloscope(harmony_embedding, cell_sample_ids,
                            dens = "GMM", dist_mat = "KL",
                            r = 10000, 
                            BPPARAM= BiocParallel::SerialParam(RNGseed=2))

dist_mat_KNN_har_sample = gloscope(harmony_embedding, cell_sample_ids,
                            dens = "KNN", dist_mat = "KL",
                            r = 10000, k=50,
                            BPPARAM= BiocParallel::SerialParam(RNGseed=2))


load("../../results/BatchStudy/rerun_harmonyallsample.Rda"))
harmony_site_embedding <- plot_df[,paste0("harmony_",1:10)]
dist_mat_GMM_har_site = gloscope(harmony_site_embedding, cell_sample_ids,
                                   dens = "GMM", dist_mat = "KL",
                                   r = 10000, 
                                   BPPARAM= BiocParallel::SerialParam(RNGseed=2))

dist_mat_KNN_har_site = gloscope(harmony_site_embedding, cell_sample_ids,
                                   dens = "KNN", dist_mat = "KL",
                                   r = 10000, k=50,
                                   BPPARAM= BiocParallel::SerialParam(RNGseed=2))

save(dist_mat_GMM_PCA,dist_mat_KNN_PCA,
     dist_mat_GMM_har_sample,dist_mat_KNN_har_sample,
     dist_mat_GMM_har_site,dist_mat_KNN_har_site,
     file = "../../results/BatchStudy/distmat_sampleid_batchevalue.Rda")

