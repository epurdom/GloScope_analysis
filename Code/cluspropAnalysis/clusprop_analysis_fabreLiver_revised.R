here::i_am("Code/cluspropAnalysis/clusprop_analysis_fabreLiver_revised.R")

source(here::here("Code","cluspropAnalysis","clusprop_helpers.R"))
dataset_name <- "fabreLiver"

# Create a df with only cell-level clusterings and cell-type labels
clustering_results <- readRDS(here::here("results","Revision","cluspropFinal",
    paste0(dataset_name,"_clustered_revised.RDS")))
cluster_columns_leiden <- colnames(clustering_results)[grep("leiden_*",colnames(clustering_results),value = FALSE)]
cluster_columns_louvain <- colnames(clustering_results)[grep("louvain_*",colnames(clustering_results),value = FALSE)]
cluster_columns <- c(cluster_columns_leiden, cluster_columns_louvain)
cluster_df <- clustering_results[,cluster_columns]
# Add "True" annotations
annotation_col <- "cell_type__ontology_label"
cluster_df$annotation <- clustering_results[,annotation_col]

# Create a df with the per-sample metadata with a sorted row order
sorted_samples <- unique(as.character(sort(clustering_results$sample)))
unsorted_sample_data <- unique(clustering_results[,c("sample","patient","group","batch")])
rownames(unsorted_sample_data) <- unsorted_sample_data$sample
sorted_sample_data <- unsorted_sample_data[sorted_samples,]

options(warn=-1) # Surpress warnings
results_list_batch_orig <- list()
results_list_phenotype_orig <- list()
results_list_batch_root <- list()
results_list_phenotype_root <- list()

for (col_index in 1:ncol(cluster_df)){
    
    col_name <- colnames(cluster_df)[col_index]
    cluster_assignments <- as.factor(cluster_df[,col_index])
    num_clusters <- length(unique(cluster_assignments))
    
    method_value <- "GloScopeProp"
    if (col_name == "annotation"){
        algorithm_value <- "Annotation"
        k_value <- NA
        resolution_value <- NA
        seed_value <- 1
    } else{
        col_name_parse <- strsplit(col_name,"_")[[1]]
        algorithm_value <- col_name_parse[1]
        k_value <- sub('.', '', strsplit(col_name,"_")[[1]][2])
        resolution_value <- sub('.', '', strsplit(col_name,"_")[[1]][3])
        seed_value <- sub('.', '', strsplit(col_name,"_")[[1]][4])
    }
    
    # Standardize the GS matrix and take a root transform
    gloscope_proportion_mat <- gloscope_proportion(clustering_results$sample,
        cluster_assignments, ep = 0.5, dist_mat = "KL")
    gloscope_proportion_mat <- gloscope_proportion_mat[sorted_samples,sorted_samples] # IMPORTANT!
    gloscope_proportion_mat <- apply(gloscope_proportion_mat,c(1,2),function(x){max(x,1e-16)})
    gloscope_proportion_mat_root <- apply(gloscope_proportion_mat,c(1,2),function(x){sqrt(x)})
    diag(gloscope_proportion_mat) <- 0
    diag(gloscope_proportion_mat_root) <- 0
    
    # Get the separation statistics
    stats_batch_orig <- get_statistics(gloscope_proportion_mat, sorted_sample_data$batch)
    stats_phenotype_orig <- get_statistics(gloscope_proportion_mat, sorted_sample_data$group)
    stats_batch_root <- get_statistics(gloscope_proportion_mat_root, sorted_sample_data$batch)
    stats_phenotype_root <- get_statistics(gloscope_proportion_mat_root, sorted_sample_data$group)
    
    # Save results
    results_vector_batch_orig <- c(method = method_value, algorithm = algorithm_value,
        k = k_value, resolution = resolution_value, replicate = seed_value,
        R = stats_batch_orig$R, O2 = stats_batch_orig$O2, silhouette = stats_batch_orig$silhouette,
        num_clusters = num_clusters, separation_var = "batch", gs_transform = "none")
    results_list_batch_orig[[col_index]] <- results_vector_batch_orig
    
    results_vector_phenotype_orig <- c(method = method_value, algorithm = algorithm_value,
        k = k_value, resolution = resolution_value, replicate = seed_value,
        R = stats_phenotype_orig$R, O2 = stats_phenotype_orig$O2, silhouette = stats_phenotype_orig$silhouette,
        num_clusters = num_clusters, separation_var = "group", gs_transform = "none")
    results_list_phenotype_orig[[col_index]] <- results_vector_phenotype_orig
    
    results_vector_batch_root <- c(method = method_value, algorithm = algorithm_value,
        k = k_value, resolution = resolution_value, replicate = seed_value,
        R = stats_batch_root$R, O2 = stats_batch_root$O2, silhouette = stats_batch_root$silhouette,
        num_clusters = num_clusters, separation_var = "batch", gs_transform = "root")
    results_list_batch_root[[col_index]] <- results_vector_batch_root
    
    results_vector_phenotype_root <- c(method = method_value, algorithm = algorithm_value,
        k = k_value, resolution = resolution_value, replicate = seed_value,
        R = stats_phenotype_root$R, O2 = stats_phenotype_root$O2, silhouette = stats_phenotype_root$silhouette,
        num_clusters = num_clusters, separation_var = "group", gs_transform = "root")
    results_list_phenotype_root[[col_index]] <- results_vector_phenotype_root
}
options(warn=0)

results_df_gloscopeProp <- do.call(rbind,
    c(results_list_batch_orig, results_list_phenotype_orig, results_list_batch_root, results_list_phenotype_root))

load(here::here("results","BatchStudy","dist_mat_fabre_liver_update.Rda"))
dist_mat_GMM_PCA_orig <- dist_mat_GMM_PCA[sorted_samples, sorted_samples]
dist_mat_GMM_PCA_orig <- apply(dist_mat_GMM_PCA_orig,c(1,2),function(x){max(x,1e-16)})
dist_mat_GMM_PCA_root <- apply(dist_mat_GMM_PCA_orig,c(1,2),function(x){sqrt(x)})
diag(dist_mat_GMM_PCA_orig) <- 0
diag(dist_mat_GMM_PCA_root) <- 0

dist_mat_KNN_PCA_orig <- dist_mat_KNN_PCA[sorted_samples, sorted_samples]
dist_mat_KNN_PCA_orig <- apply(dist_mat_KNN_PCA_orig,c(1,2),function(x){max(x,1e-16)})
dist_mat_KNN_PCA_root <- apply(dist_mat_KNN_PCA_orig,c(1,2),function(x){sqrt(x)})
diag(dist_mat_KNN_PCA_orig) <- 0
diag(dist_mat_KNN_PCA_root) <- 0

# GMM + Original
full_gloscope_results_gmm_batch_orig <- get_statistics(dist_mat_GMM_PCA_orig, sorted_sample_data$batch)
results_vector_gmm_batch_orig <- c(method = "GloScopeGMM", algorithm = NA,
        k = NA, resolution = NA, replicate = 1,
        R = full_gloscope_results_gmm_batch_orig$R,
        O2 = full_gloscope_results_gmm_batch_orig$O2,
        silhouette = full_gloscope_results_gmm_batch_orig$silhouette,
        num_clusters = NA, separation_var = "batch", gs_transform = "none")

full_gloscope_results_gmm_phenotype_orig <- get_statistics(dist_mat_GMM_PCA_orig, sorted_sample_data$group)
results_vector_gmm_phenotype_orig <- c(method = "GloScopeGMM", algorithm = NA,
        k = NA, resolution = NA, replicate = 1,
        R = full_gloscope_results_gmm_phenotype_orig$R,
        O2 = full_gloscope_results_gmm_phenotype_orig$O2,
        silhouette = full_gloscope_results_gmm_phenotype_orig$silhouette,
        num_clusters = NA, separation_var = "group", gs_transform = "none")

# GMM + Root
full_gloscope_results_gmm_batch_root <- get_statistics(dist_mat_GMM_PCA_root, sorted_sample_data$batch)
results_vector_gmm_batch_root <- c(method = "GloScopeGMM", algorithm = NA,
        k = NA, resolution = NA, replicate = 1,
        R = full_gloscope_results_gmm_batch_root$R,
        O2 = full_gloscope_results_gmm_batch_root$O2,
        silhouette = full_gloscope_results_gmm_batch_root$silhouette,
        num_clusters = NA, separation_var = "batch", gs_transform = "root")

full_gloscope_results_gmm_phenotype_root <- get_statistics(dist_mat_GMM_PCA_root, sorted_sample_data$group)
results_vector_gmm_phenotype_root <- c(method = "GloScopeGMM", algorithm = NA,
        k = NA, resolution = NA, replicate = 1,
        R = full_gloscope_results_gmm_phenotype_root$R,
        O2 = full_gloscope_results_gmm_phenotype_root$O2,
        silhouette = full_gloscope_results_gmm_phenotype_root$silhouette,
        num_clusters = NA, separation_var = "group", gs_transform = "root")

# KNN + Original
full_gloscope_results_knn_batch_orig <- get_statistics(dist_mat_KNN_PCA_orig, sorted_sample_data$batch)
results_vector_knn_batch_orig <- c(method = "GloScopeKNN", algorithm = NA,
        k = NA, resolution = NA, replicate = 1,
        R = full_gloscope_results_knn_batch_orig$R,
        O2 = full_gloscope_results_knn_batch_orig$O2,
        silhouette = full_gloscope_results_knn_batch_orig$silhouette,
        num_clusters = NA, separation_var = "batch", gs_transform = "none")

full_gloscope_results_knn_phenotype_orig <- get_statistics(dist_mat_KNN_PCA_orig, sorted_sample_data$group)
results_vector_knn_phenotype_orig <- c(method = "GloScopeKNN", algorithm = NA,
        k = NA, resolution = NA, replicate = 1,
        R = full_gloscope_results_knn_phenotype_orig$R,
        O2 = full_gloscope_results_knn_phenotype_orig$O2,
        silhouette = full_gloscope_results_knn_phenotype_orig$silhouette,
        num_clusters = NA, separation_var = "group", gs_transform = "none")

# KNN + Root
full_gloscope_results_knn_batch_root <- get_statistics(dist_mat_KNN_PCA_root, sorted_sample_data$batch)
results_vector_knn_batch_root <- c(method = "GloScopeKNN", algorithm = NA,
        k = NA, resolution = NA, replicate = 1,
        R = full_gloscope_results_knn_batch_root$R,
        O2 = full_gloscope_results_knn_batch_root$O2,
        silhouette = full_gloscope_results_knn_batch_root$silhouette,
        num_clusters = NA, separation_var = "batch", gs_transform = "root")

full_gloscope_results_knn_phenotype_root <- get_statistics(dist_mat_KNN_PCA_root, sorted_sample_data$group)
results_vector_knn_phenotype_root <- c(method = "GloScopeKNN", algorithm = NA,
        k = NA, resolution = NA, replicate = 1,
        R = full_gloscope_results_knn_phenotype_root$R,
        O2 = full_gloscope_results_knn_phenotype_root$O2,
        silhouette = full_gloscope_results_knn_phenotype_root$silhouette,
        num_clusters = NA, separation_var = "group", gs_transform = "root")

results_df_gloscopeFull <- rbind(results_vector_gmm_batch_orig, results_vector_gmm_phenotype_orig,
    results_vector_gmm_batch_root, results_vector_gmm_phenotype_root,
    results_vector_knn_batch_orig, results_vector_knn_phenotype_orig,
    results_vector_knn_batch_root, results_vector_knn_phenotype_root)
results_df_gloscopeAll <- rbind(results_df_gloscopeProp, results_df_gloscopeFull)

library(sceasy)
library(reticulate)
use_condaenv("PILOT")
pl <- import("PILOT")
sc <- import("scanpy")

if(!file.exists(here::here("results","Revision","cluspropFinal",paste0(dataset_name,".h5ad")))){
    seurat_object <- readRDS(here::here("data","Processed_Datasets",dataset_name,
        paste0(dataset_name,"_default"), paste0(dataset_name,"_default.Rds")))
    seurat_object[[]] <- clustering_results
    seurat_object[["RNA"]] <- as(seurat_object[["RNA"]], "Assay")

    sceasy::convertFormat(seurat_object, from="seurat", to="anndata",
        outFile=here::here("results","Revision","cluspropFinal",paste0(dataset_name,".h5ad")))
}
anndata <- sc$read_h5ad(here::here("results","Revision","cluspropFinal",paste0(dataset_name,".h5ad")))

results_list_batch <- list()
results_list_phenotype <- list()

for (col_index in 1:ncol(cluster_df)){
    
    col_name <- colnames(cluster_df)[col_index]
    cluster_assignments <- as.factor(cluster_df[,col_index])
    num_clusters <- length(unique(cluster_assignments))
    
    anndata$obs[col_name] <- cluster_assignments
    
    method_value <- "PILOT"
    if (col_name == "annotation"){
        algorithm_value <- "Annotation"
        k_value <- NA
        resolution_value <- NA
        seed_value <- 1
    } else{
        col_name_parse <- strsplit(col_name,"_")[[1]]
        algorithm_value <- col_name_parse[1]
        k_value <- sub('.', '', strsplit(col_name,"_")[[1]][2])
        resolution_value <- sub('.', '', strsplit(col_name,"_")[[1]][3])
        seed_value <- sub('.', '', strsplit(col_name,"_")[[1]][4])
    }
        
    pl$tl$wasserstein_distance(
        anndata,
        emb_matrix = 'X_pca',
        clusters_col = col_name,
        sample_col = 'sample',
        status = 'group',
        )
    
    pilot_mat <- anndata$uns$EMD_df
    pilot_mat <- pilot_mat[sorted_samples,sorted_samples]
    pilot_mat <- as.matrix(pilot_mat)
    
    pilot_stats_batch <- get_statistics(pilot_mat, sorted_sample_data$batch)
    results_vector_batch <- c(method = method_value, algorithm = algorithm_value,
        k = k_value, resolution = resolution_value, replicate = seed_value,
        R = pilot_stats_batch$R, O2 = pilot_stats_batch$O2, silhouette = pilot_stats_batch$silhouette,
        num_clusters = num_clusters, separation_var = "batch", gs_transform = NA)
    
    pilot_stats_phenotype <- get_statistics(pilot_mat, sorted_sample_data$group)
    results_vector_phenotype <- c(method = method_value, algorithm = algorithm_value,
        k = k_value, resolution = resolution_value, replicate = seed_value,
        R = pilot_stats_phenotype$R, O2 = pilot_stats_phenotype$O2, silhouette = pilot_stats_phenotype$silhouette,
        num_clusters = num_clusters, separation_var = "group", gs_transform = NA)
    
    results_list_batch[[col_index]] <- results_vector_batch
    results_list_phenotype[[col_index]] <- results_vector_phenotype
}

results_df_PILOT <- do.call(rbind,c(results_list_batch,results_list_phenotype))
results_df_all <- rbind(results_df_gloscopeAll, results_df_PILOT)
saveRDS(results_df_all, here::here("results","Revision","cluspropFinal",paste0(dataset_name,"_results.RDS")))