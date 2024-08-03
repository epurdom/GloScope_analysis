here::i_am("Code/Revision/timing/time_pca_perezLupus")

library(GloScope)
library(SingleCellExperiment)                                                
library(scuttle) 
library(scater)
library(HDF5Array)
library(umap)
library(igraph)
library(scran)

dataset_name <- "perezLupus"
program_name <- "default"
sample_covar <- "sample"

set.seed(2)
timing_vector <- c()

# Covert raw files into a Seurat object

to_seurat_start <- proc.time()[3]

#system(paste0("python3.10 ",here::here("Code","Process_Datasets","Programs","perezLupus_default.py")))
raw_file_path <- here::here("data","Processed_Datasets","perezLupus","raw_files","scRNA_raw_cleaned.h5ad")
log_object <- HDF5Array::H5ADMatrix(raw_file_path, "data")
count_object <- HDF5Array::H5ADMatrix(raw_file_path, "counts")
   
metadata_df <- readRDS(here::here("data", "Processed_Datasets", "perezLupus", "metadata", "perezLupus_default_metadata.Rds"))

sce_object <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = count_object, log = log_object), colData = metadata_df)

sce_object@assays@data@listData$counts@seed@filepath <- "/scratch/users/singlecell/pop_data/data/Processed_Datasets/perezLupus/raw_files/scRNA_raw_cleaned.h5ad"
sce_object@assays@data@listData$log@seed@filepath <- "/scratch/users/singlecell/pop_data/data/Processed_Datasets/perezLupus/raw_files/scRNA_raw_cleaned.h5ad"
#sce_object@assays@data@listData$logcounts@seed@seed@seed@seed@seed@filepath <- "/scratch/users/singlecell/pop_data/data/Processed_Datasets/perezLupus/raw_files/scRNA_raw_cleaned.h5ad"

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
