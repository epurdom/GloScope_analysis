here::i_am("Code/timingAnalysis/yaoMouseBrain_timing_simple.R")

library(GloScope)
library(SingleCellExperiment)                                                
library(scuttle) 
library(scater)
library(HDF5Array)
library(umap)
library(igraph)
library(scran)
remotes::install_github("drighelli/AllenInstituteBrainData")

dataset_name <- "yaoMouseBrain"
program_name <- "default"
sample_covar <- "sample"

set.seed(2)
timing_vector <- c()

# Covert raw files into a Seurat object

to_seurat_start <- proc.time()[3]

metadata_df <- readRDS(here::here("data", "Processed_Datasets", "yaoMouseBrain", "metadata", "yaoMouseBrain_default_metadata.Rds"))

AllenInstituteBrainData:::.initCache(here::here("data","Processed_Datasets","yaoMouseBrain","sce_files"))
sce_object <- AllenInstituteBrainData::AllenInstituteBrainData("Allen_Mouse_2020")
sce_object$sample <- interaction(sce_object$donor_label,sce_object$region_label,sep="_")

to_seurat_end <- proc.time()[3]

to_seurat_time <- to_seurat_end - to_seurat_start
timing_vector["to_seurat"] <- to_seurat_time

# Normalize data

normalize_start <- proc.time()[3]

sce_object <- scuttle::logNormCounts(sce_object)

normalize_end <- proc.time()[3]

normalize_time <- normalize_end - normalize_start
timing_vector["normalize"] <- normalize_time

# Filter to HVGS

hvg_start <- proc.time()[3]

variable_genes <- scran::modelGeneVar(sce_object)                 
selected_genes <- scran::getTopHVGs(variable_genes, n = 2000)     

hvg_end <- proc.time()[3]
hvg_time <- hvg_end - hvg_start
timing_vector["filter"] <- hvg_time

# Run PCA

pca_start <- proc.time()[3]

sce_object <- scater::runPCA(sce_object,scale = TRUE, subset_row = selected_genes)

pca_end <- proc.time()[3]

pca_time <- pca_end - pca_start
timing_vector["pca"] <- pca_time

save_path <- paste0(here::here("Code","timingAnalysis","timing_results",paste0(dataset_name,"_pca.RDS")))
saveRDS(timing_vector,save_path)
