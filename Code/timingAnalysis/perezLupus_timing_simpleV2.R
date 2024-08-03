here::i_am("Code/timingAnalysis/perezLupus_timing_simple.R")

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

timing_vector["filter"] <- NA

# Scale and Center

timing_vector["scale"] <- NA

# Run PCA

pca_start <- proc.time()[3]

sce_object <- runPCA(sce_object)

pca_end <- proc.time()[3]

pca_time <- pca_end - pca_start
timing_vector["pca"] <- pca_time

# Run GloScope

embedding_df <- sce_object@int_colData@listData$reducedDims$PCA[,1:10]
sample_ids <- sce_object$sample
print(BiocParallel::SerialParam())

gloscope_gmm_start <- proc.time()[3]

gloscope_gmm  <- GloScope::gloscope(embedding_df, sample_ids, dens = "GMM",
    BPPARAM = BiocParallel::SerialParam())

gloscope_gmm_end <- proc.time()[3]
gloscope_gmm_time <- gloscope_gmm_end - gloscope_gmm_start
timing_vector["gloscope_gmm"] <- gloscope_gmm_time

gloscope_knn_start <- proc.time()[3]

gloscope_knn  <- GloScope::gloscope(embedding_df, sample_ids, dens = "KNN",
    BPPARAM = BiocParallel::SerialParam())

gloscope_knn_end <- proc.time()[3]
gloscope_knn_time <- gloscope_knn_end - gloscope_knn_start
timing_vector["gloscope_knn"] <- gloscope_knn_time

# Run UMAP

umap_start <- proc.time()[3]

pca_umap <- umap(embedding_df[,1:10])

umap_end <- proc.time()[3]

umap_time <- umap_end - umap_start
timing_vector["umap"] <- umap_time

# Compute NN graph

nn_start <- proc.time()[3]

k_values <- c(5, 20, 100)
reducedDim(sce_object, "PCAsubset") <- reducedDim(sce_object, "PCA")[,1:10]

graph_k20 <- scran::buildKNNGraph(sce_object, use.dimred="PCAsubset", k = 20)

nn_end <- proc.time()[3]

nn_time <- nn_end - nn_start
timing_vector["nn"] <- nn_time

# Find Clusters

cluster_start <- proc.time()[3]

modularity_values <- c(0.1, 0.8, 2)


membership_assignments <- igraph::cluster_louvain(graph_k20,
                resolution = 1)$membership


cluster_end <- proc.time()[3]

cluster_time <- cluster_end - cluster_start
timing_vector["cluster"] <- cluster_time

num_cells <- ncol(sce_object)
num_genes <- nrow(sce_object)
num_samples <- length(unique(sce_object$sample))

timing_vector["num_genes"] <- num_genes
timing_vector["num_cells"] <- num_cells
timing_vector["num_samples"] <- num_samples


save_path <- here::here("Code","timingAnalysis","timing_results",paste0(dataset_name,"_simple.RDS"))
saveRDS(timing_vector,save_path)
