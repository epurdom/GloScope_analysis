here::i_am("Code/cluspropAnalysis/clusprop_clustering_fabreLung_revised.R")

library(igraph)
library(Seurat)

dataset_name <- "fabreLung"
seurat_object <- readRDS(here::here("data","Processed_Datasets",
    dataset_name,paste0(dataset_name,"_default"),
    paste0(dataset_name,"_default.Rds")))

seurat_object <- Seurat::FindNeighbors(seurat_object, k.param = 5, graph.name = "knn_graph_k5")
seurat_object <- Seurat::FindNeighbors(seurat_object, k.param = 20, graph.name = "knn_graph_k20")
seurat_object <- Seurat::FindNeighbors(seurat_object, k.param = 100, graph.name = "knn_graph_k100")

igraph_k5 <- igraph::graph_from_adjacency_matrix(seurat_object@graphs$knn_graph_k5,
    mode = "undirected")
igraph_k20 <- igraph::graph_from_adjacency_matrix(seurat_object@graphs$knn_graph_k20,
    mode = "undirected")
igraph_k100 <- igraph::graph_from_adjacency_matrix(seurat_object@graphs$knn_graph_k100,
    mode = "undirected")

k_vec <- c(5, 20, 100)
resolution_vec_leiden <- c(1e-4, 5e-5, 1e-5, 5e-6)
resolution_vec_louvain <- c(0.1, 0.8, 2, 5)
num_replicates <- 20

save_path <- paste0("Code/cluspropAnalysis/clustered_data/",dataset_name,"_clustered_revised.RDS")
for (k_value in k_vec){
    print(k_value)
    graph_name <- paste0("igraph_k",k_value)
    for (replicate_index in 1:num_replicates){
        print(replicate_index)
        for (resolution_value_leiden in resolution_vec_leiden){
            # Run Leiden
            id_string <- paste0("leiden_k",k_value,"_r",resolution_value_leiden,"_s",replicate_index)
            set.seed(replicate_index)
            membership_assignments <- igraph::cluster_leiden(get(graph_name),
                resolution_parameter = resolution_value_leiden)$membership
            seurat_object[[]][,id_string] <- membership_assignments
        }
        for (resolution_value_louvain in resolution_vec_louvain){
            # Run Louvian
            id_string <- paste0("louvain_k",k_value,"_r",resolution_value_louvain,"_s",replicate_index)
            set.seed(replicate_index)
            membership_assignments <- igraph::cluster_louvain(get(graph_name),
                resolution = resolution_value_louvain)$membership
            seurat_object[[]][,id_string] <- membership_assignments
        }
        saveRDS(seurat_object[[]],here::here(save_path))
    }
}

saveRDS(seurat_object[[]],here::here(save_path))