here::i_am("Code/timingAnalysis/time_pca_perezLupus2.R")

library(GloScope)
library(SingleCellExperiment)                                                
library(scater)
library(HDF5Array)

dataset_name <- "perezLupus"
program_name <- "default"
sample_covar <- "sample"

set.seed(2)

sce_object <- readRDS(here::here("data","Processed_Datasets","perezLupus",
    "perezLupus_default", "perezLupus_default.Rds"))

gene_list <- rownames(attributes(reducedDim(sce_object, "PCA"))$`rotation`)

sce_object@assays@data@listData$counts@seed@filepath <- "/scratch/users/singlecell/pop_data/data/Processed_Datasets/perezLupus/raw_files/scRNA_raw_cleaned.h5ad"
sce_object@assays@data@listData$log@seed@filepath <- "/scratch/users/singlecell/pop_data/data/Processed_Datasets/perezLupus/raw_files/scRNA_raw_cleaned.h5ad"
sce_object@assays@data@listData$logcounts@seed@seed@seed@seed@seed@filepath <- "/scratch/users/singlecell/pop_data/data/Processed_Datasets/perezLupus/raw_files/scRNA_raw_cleaned.h5ad"

pca_start <- proc.time()[3]

sce_object <- scater::runPCA(sce_object,scale = TRUE, subset_row = gene_list)

pca_end <- proc.time()[3]

pca_time <- pca_end - pca_start
print(pca_time)