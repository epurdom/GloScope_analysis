here::i_am("Code/timingAnalysis/stephensonCOVIDPBMC_timing_simple.R")

library(GloScope)
library(Seurat)

dataset_name <- "stephensonCOVIDPBMC"
program_name <- "default"
sample_covar <- "sample_id"

set.seed(2)
timing_vector <- c()

# Covert raw files into a Seurat object

to_seurat_start <- proc.time()[3]

seurat_object <- SeuratDisk::LoadH5Seurat(here::here("data","Processed_Datasets","stephensonCOVIDPBMC","raw_files","covid_portal_210320_with_raw.h5seurat"))


to_seurat_end <- proc.time()[3]

to_seurat_time <- to_seurat_end - to_seurat_start
timing_vector["to_seurat"] <- to_seurat_time

# Normalize data

normalize_start <- proc.time()[3]

seurat_object <- Seurat::NormalizeData(seurat_object)

normalize_end <- proc.time()[3]

normalize_time <- normalize_end - normalize_start
timing_vector["normalize"] <- normalize_time

# Filter to HVGS

hvg_start <- proc.time()[3]

seurat_object <- Seurat::FindVariableFeatures(seurat_object)

hvg_end <- proc.time()[3]

hvg_time <- hvg_end - hvg_start
timing_vector["filter"] <- hvg_time

# Scale and Center

scale_start <- proc.time()[3]

seurat_object <- Seurat::ScaleData(seurat_object)

scale_end <- proc.time()[3]

scale_time <- scale_end - scale_start
timing_vector["scale"] <- scale_time

# Run PCA

pca_start <- proc.time()[3]

seurat_object <- Seurat::RunPCA(seurat_object)

pca_end <- proc.time()[3]

pca_time <- pca_end - pca_start
timing_vector["pca"] <- pca_time

save_path <- paste0(here::here("Code","timingAnalysis","timing_results",paste0(dataset_name,"_pca.RDS")))
saveRDS(timing_vector,save_path)