devtools::load_all(("~/popPackage"))
library(dplyr)
library(MASS)
library(SingleCellExperiment)
sce <- readRDS("../../data/Processed_Datasets/perezLupus/perezLupus_default/perezLupus_default.Rds")
#load("../../results/Perez2022/batchid_harmony.Rda")
cell_sample_ids <- sce$sample
pca_embedding <- reducedDim(sce)
harmony_embedding <-reducedDim(sce, "HARMONY")
load("../../results/Perez2022/harmonySam_pipeline.Rda")


dist_mat_GMM_PCA = gloscope(pca_embedding[,1:10], cell_sample_ids,
                            dens = "GMM", dist_mat = "KL",
                            r = 10000, 
                            BPPARAM= BiocParallel::MulticoreParam(2,RNGseed = 1))



dist_mat_GMM_har_sample = gloscope(harmonized_pcs[,1:10], cell_sample_ids,
                                   dens = "GMM", dist_mat = "KL",
                                   r = 10000, 
                                   BPPARAM= BiocParallel::SerialParam(RNGseed=2))


dist_mat_GMM_har_batch = gloscope(harmony_embedding[,1:10], cell_sample_ids,
                                  dens = "GMM", dist_mat = "KL",
                                  r = 10000, 
                                  BPPARAM= BiocParallel::MulticoreParam(2,RNGseed = 1))



dist_mat_KNN_PCA = gloscope(pca_embedding[,1:10], cell_sample_ids,
                            dens = "KNN", dist_mat = "KL",
                            r = 10000, 
                            BPPARAM= BiocParallel::MulticoreParam(2,RNGseed = 1))


dist_mat_KNN_har_batch = gloscope(harmony_embedding[,1:10], cell_sample_ids,
                                  dens = "KNN", dist_mat = "KL",
                                  r = 10000, 
                                  BPPARAM= BiocParallel::MulticoreParam(2,RNGseed = 1))

dist_mat_KNN_har_sample = gloscope(harmonized_pcs[,1:10], cell_sample_ids,
                                   dens = "KNN", dist_mat = "KL",
                                   r = 10000, 
                                   BPPARAM= BiocParallel::MulticoreParam(2,RNGseed = 1))


save(dist_mat_GMM_PCA,dist_mat_GMM_har_sample, dist_mat_GMM_har_batch,
     dist_mat_KNN_PCA,dist_mat_KNN_har_sample,dist_mat_KNN_har_batch,
     file = "../../results/BatchStudy/dist_mat_Perez2022.Rda")
