devtools::load_all(here::here("..","popPackage"))
library(dplyr)
library(MASS)
library(mclust)
latent <- read.csv("../../results/COVID_143/scvi_nobatch.csv")[,-1]
latent_batch <- read.csv("../../results/COVID_143/scvi_test.csv")[,-1]
latent_batch_pat <- read.csv("../../results/COVID_143/scvi_sam.csv")[,-1]
meta <- readRDS("../../data/Processed_Datasets/stephensonCOVIDPBMC/metadata/stephensonCOVIDPBMC_default_metadata.Rds")


colnames(latent_batch) <- paste0("scsite_", 1:ncol(latent_batch))
colnames(latent) <- paste0("scvi_", 1:ncol(latent))
colnames(latent_batch_pat) <- paste0("scsamp_", 1:ncol(latent_batch))

cell_sample_ids <- plot_df$sample


dist_mat_GMM_scvi =  gloscope(latent, cell_sample_ids,
                              dens = "GMM", dist_mat = "KL",
                              r = 10000, 
                              BPPARAM= BiocParallel::SerialParam(RNGseed=2))

dist_mat_KNN_scvi = gloscope(latent, cell_sample_ids,
                             dens = "KNN", dist_mat = "KL",
                             r = 10000, 
                             BPPARAM= BiocParallel::SerialParam(RNGseed=2))

dist_mat_GMM_scvi_site = gloscope(latent_batch, cell_sample_ids,
                                  dens = "GMM", dist_mat = "KL",
                                  r = 10000, 
                                  BPPARAM= BiocParallel::SerialParam(RNGseed=2))

dist_mat_KNN_scvi_site= gloscope(latent_batch, cell_sample_ids,
                                 dens = "KNN", dist_mat = "KL",
                                 r = 10000, 
                                 BPPARAM= BiocParallel::SerialParam(RNGseed=2))


dist_mat_GMM_scvi_samp = gloscope(latent_batch_pat, cell_sample_ids,
                                  dens = "GMM", dist_mat = "KL",
                                  r = 10000, 
                                  BPPARAM= BiocParallel::SerialParam(RNGseed=2))


dist_mat_KNN_scvi_samp = gloscope(latent_batch_pat, cell_sample_ids,
                                  dens = "KNN", dist_mat = "KL",
                                  r = 10000, 
                                  BPPARAM= BiocParallel::SerialParam(RNGseed=2))

save(
  dist_mat_GMM_scvi, dist_mat_KNN_scvi,
     dist_mat_GMM_scvi_site, dist_mat_KNN_scvi_site,
     dist_mat_GMM_scvi_samp, dist_mat_KNN_scvi_samp,
     file = "../../results/BatchStudy/scvi_distmat_covid143.Rda")


