here::i_am("Code/timingAnalysis/chengRash12_timingSerial.R")

library(GloScope)
library(Seurat)

dataset_name <- "chengRash12"
program_name <- "default"
sample_covar <- "sample"

set.seed(2)
timing_vector <- c()

# Covert raw files into a Seurat object

to_seurat_start <- proc.time()[3]

count_data <- data.table::fread(here::here("data","Processed_Datasets",dataset_name,"raw_files","EGAF00002195824","filtered_counts.csv.gz"))
count_data <- tibble::column_to_rownames(count_data,var="V1")

metadata_df <- read.table(here::here("data","Processed_Datasets",dataset_name,"raw_files","EGAF00002195821","publication_celldata.csv"),header=T,row.names=1,sep=",") 

seurat_object <- Seurat::CreateSeuratObject(counts=count_data,project="timing_test",assay="RNA",meta.data=metadata_df)

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

# Run GloScope

embedding_df <- seurat_object@reductions$pca@cell.embeddings
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

k_values <- c(5, 20, 100)

for (k in k_values){
    graph_name <- paste0("graph_",k)
    seurat_object <- Seurat::FindNeighbors(seurat_object, k.param = k, graph.name = graph_name)
}

nn_end <- proc.time()[3]

nn_time <- nn_end - nn_start
timing_vector["nn"] <- nn_time

# Find Clusters

cluster_start <- proc.time()[3]

modularity_values <- c(0.1, 0.8, 2)

for (k in k_values){
    graph_name <- paste0("graph_",k)
    for (m in modularity_values){
        cluster_name <- paste0("graph_",k,"_modularity_",m)
        seurat_object <- Seurat::FindClusters(seurat_object, algorithm = 1, resolution = m,
            graph.name = graph_name, cluster.name = cluster_name)
    }
}


cluster_end <- proc.time()[3]

cluster_time <- cluster_end - cluster_start
timing_vector["cluster"] <- cluster_time

num_cells <- ncol(seurat_object)
num_genes <- nrow(seurat_object)
num_samples <- length(unique(seurat_object[[]][,sample_covar]))

timing_vector["num_genes"] <- num_genes
timing_vector["num_cells"] <- num_cells
timing_vector["num_samples"] <- num_samples

save_path <- paste0("Code/timingAnalysis/timing_results/",dataset_name,"_serial.RDS")
saveRDS(timing_vector,save_path)