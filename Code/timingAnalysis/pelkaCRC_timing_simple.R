here::i_am("Code/timingAnalysis/pelkaCRC_timing_simple.R")

library(GloScope)
library(Seurat)

dataset_name <- "pelkaCRC"
program_name <- "default"
sample_covar <- "PatientTypeID"

set.seed(2)
timing_vector <- c()

# Covert raw files into a Seurat object

to_seurat_start <- proc.time()[3]

count_data <- Seurat::Read10X_h5(here::here("data","Processed_Datasets","pelkaCRC","raw_files","GSE178341_crc10x_full_c295v4_submit.h5"))

  
metadata_load <- read.csv(here::here("data","Processed_Datasets","pelkaCRC","raw_files","meta.tsv"),sep="\t",row.names=1)[-1,]
cell_types_df <- read.csv(here::here("data","Processed_Datasets","pelkaCRC","raw_files","GSE178341_crc10x_full_c295v4_submit_cluster.csv"),row.names=1)
metadata_df <- merge(cell_types_df, metadata_load, by = "row.names")
rownames(metadata_df) <- metadata_df[,1]
metadata_df[,1] <- NULL # remove now duplicate column in-place

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
low_cell_samples <-  c("C119_N") #c("C114_N","C115_N","C119_N","C133_N","C135_N","C151_T")
low_cell_mask <- !(seurat_object[[]]$PatientTypeID %in% low_cell_samples)
filtered_names <- colnames(seurat_object)[low_cell_mask]
seurat_object <- seurat_object[,colnames(seurat_object) %in% filtered_names]

	# Two samples have 10X V2 and V3 sequencing. We only keep the V3 data
sample_mask <- (seurat_object[[]]$donor_id != "C162")
tech_mask <- (seurat_object[[]]$library_preparation_protocol__ontology_label != "10X 3' v2 sequencing")
filter_mask <- sample_mask | tech_mask
filtered_names <- colnames(seurat_object)[filter_mask]
seurat_object <- seurat_object[,colnames(seurat_object) %in% filtered_names]

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
