here::i_am("Code/cluspropAnalysis/clusprop_clustering_perezLupus_revised.R")

library(igraph)
library(scran)
library(SingleCellExperiment)

dataset_name <- "perezLupus"

sce <- readRDS(here::here("data","Processed_Datasets","perezLupus","perezLupus_default","perezLupus_default.Rds"))
reducedDim(sce, "PCAsubset") <- reducedDim(sce, "PCA")[,1:10] # Match Seurat

graph_k5 <- scran::buildKNNGraph(sce, use.dimred="PCAsubset", k = 5)
graph_k20 <- scran::buildKNNGraph(sce, use.dimred="PCAsubset", k = 20)
graph_k100 <- scran::buildKNNGraph(sce, use.dimred="PCAsubset", k = 100)

k_vec <- c(5,20,100)
resolution_vec_leiden <- c(1e-4, 5e-5, 1e-5, 5e-6)
resolution_vec_louvain <- c(0.1, 0.8, 2, 5)
num_replicates <- 20

save_path <- paste0("Code/cluspropAnalysis/clustered_data/",dataset_name,"_clustered_revised.RDS")
for (k_value in k_vec){
    print(k_value)
    graph_name <- paste0("graph_k",k_value)
    for (replicate_index in 1:num_replicates){
        print(replicate_index)
        for (resolution_value_leiden in resolution_vec_leiden){
            # Run Leiden
            id_string <- paste0("leiden_k",k_value,"_r",resolution_value_leiden,"_s",replicate_index)
            set.seed(replicate_index)
            membership_assignments <- igraph::cluster_leiden(get(graph_name),
                resolution_parameter = resolution_value_leiden)$membership
            colData(sce)[[id_string]] <- membership_assignments
        }
        for (resolution_value_louvain in resolution_vec_louvain){
            # Run Louvian
            id_string <- paste0("louvain_k",k_value,"_r",resolution_value_louvain,"_s",replicate_index)
            set.seed(replicate_index)
            membership_assignments <- igraph::cluster_louvain(get(graph_name),
                resolution = resolution_value_louvain)$membership
            colData(sce)[[id_string]] <- membership_assignments
        }
        saveRDS(colData(sce),here::here(save_path))
    }
}