devtools::load_all(("~/popPackage"))
library(dplyr)
library(MASS)
library(SingleCellExperiment)
sce <- readRDS("../../data/Processed_Datasets/perezLupus/perezLupus_default/perezLupus_default.Rds")
#load("../../results/Perez2022/batchid_harmony.Rda")
cell_sample_ids <- sce$sample


scvi_batch = read.csv("../../results/Perez2022/scvi_batch_update.csv")[,-1]
colnames(scvi_batch) = paste0("scsamp_", 1:ncol(scvi_batch))
scvi_sample = read.csv("../../results/Perez2022/scvi_batch_sample.csv")[,-1]
colnames(scvi_sample) = paste0("scvi_", 1:ncol(scvi_sample))
scvi_nobatch = read.csv("../../results/Perez2022/scvi_nobatch.csv")[,-1]
colnames(scvi_nobatch) = paste0("scvi_", 1:ncol(scvi_nobatch))




dist_mat_GMM_scvi = gloscope(scvi_nobatch, cell_sample_ids,
                             dens = "GMM", dist_mat = "KL",
                             r = 10000, 
                             BPPARAM= BiocParallel::MulticoreParam(2,RNGseed = 1))


dist_mat_GMM_scvi_sample = gloscope(scvi_sample, cell_sample_ids,
                                    dens = "GMM", dist_mat = "KL",
                                    r = 10000, 
                                    BPPARAM= BiocParallel::MulticoreParam(2,RNGseed = 1))

dist_mat_GMM_scvi_batch = gloscope(scvi_batch, cell_sample_ids,
                                   dens = "GMM", dist_mat = "KL",
                                   r = 10000, 
                                   BPPARAM= BiocParallel::MulticoreParam(2,RNGseed = 1))





dist_mat_KNN_scvi_sample = gloscope(scvi_sample, cell_sample_ids,
                                    dens = "KNN", dist_mat = "KL",
                                    r = 10000, 
                                    BPPARAM= BiocParallel::MulticoreParam(2,RNGseed = 1))

dist_mat_KNN_scvi = gloscope(scvi_nobatch, cell_sample_ids,
                             dens = "KNN", dist_mat = "KL",
                             r = 10000, 
                             BPPARAM= BiocParallel::MulticoreParam(2,RNGseed = 1))



dist_mat_KNN_scvi_batch = gloscope(scvi_batch, cell_sample_ids,
                                   dens = "KNN", dist_mat = "KL",
                                   r = 10000, 
                                   BPPARAM= BiocParallel::MulticoreParam(2,RNGseed = 1))




save(dist_mat_GMM_scvi, dist_mat_GMM_scvi_batch,dist_mat_GMM_scvi_sample,
     dist_mat_KNN_scvi, dist_mat_KNN_scvi_batch,dist_mat_KNN_scvi_sample,
     file = "../../results/BatchStudy/dist_mat_Perez2022_scvi.Rda")
