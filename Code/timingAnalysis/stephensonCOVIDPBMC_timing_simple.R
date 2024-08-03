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

timing_vector["normalize"] <- NA

# Filter to HVGS

timing_vector["filter"] <- NA

# Scale and Center

timing_vector["scale"] <- NA

# Run PCA

timing_vector["pca"] <- NA

# Run GloScope

embedding_df <- seurat_object@reductions$pca@cell.embeddings[,1:10]
sample_ids <- seurat_object@meta.data[,sample_covar]
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

seurat_object <- Seurat::RunUMAP(seurat_object, reduction = "pca", dims=1:10)

umap_end <- proc.time()[3]

umap_time <- umap_end - umap_start
timing_vector["umap"] <- umap_time

# Compute NN graph

nn_start <- proc.time()[3]


    graph_name <- paste0("graph_",20)
    seurat_object <- Seurat::FindNeighbors(seurat_object, k.param = 20, graph.name = graph_name)

nn_end <- proc.time()[3]

nn_time <- nn_end - nn_start
timing_vector["nn"] <- nn_time

# Find Clusters

cluster_start <- proc.time()[3]

modularity_values <- c(0.1, 0.8, 2)


    graph_name <- paste0("graph_",20)
        cluster_name <- paste0("graph_",20,"_modularity_",0.8)
        seurat_object <- Seurat::FindClusters(seurat_object, algorithm = 1, resolution = 0.8,
            graph.name = graph_name, cluster.name = cluster_name)


cluster_end <- proc.time()[3]

cluster_time <- cluster_end - cluster_start
timing_vector["cluster"] <- cluster_time

num_cells <- ncol(seurat_object)
num_genes <- nrow(seurat_object)
num_samples <- length(unique(seurat_object[[]][,sample_covar]))

timing_vector["num_genes"] <- num_genes
timing_vector["num_cells"] <- num_cells
timing_vector["num_samples"] <- num_samples

save_path <- here::here("Code","timingAnalysis","timing_results",paste0(dataset_name,"_simple.RDS"))
saveRDS(timing_vector,save_path)