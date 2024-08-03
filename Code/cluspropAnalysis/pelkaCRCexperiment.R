here::i_am("Code/cluspropAnalysis/pelkaCRCexperiment.R")

source(here::here("Code","cluspropAnalysis","clusprop_helpers.R"))

clustered_data <- readRDS(here::here("results","Revision","cluspropFinal","pelkaCRC_clustered_revised.RDS")) 
annotation_col <- "cell_type"
clustered_data$annotation <- clustered_data[,annotation_col]

partition_sizes <- c(2,5,10)
num_random_partitions <- 5
partition_name_vec <- c(sapply(partition_sizes, function(partition_size){
    paste0("partition_",partition_size,"_",1:num_random_partitions)}))

cluster_columns_leiden <- colnames(clustered_data)[grep("leiden_*",colnames(clustered_data),value = FALSE)]
cluster_columns_louvain <- colnames(clustered_data)[grep("louvain_*",colnames(clustered_data),value = FALSE)]
cluster_columns <- c(cluster_columns_leiden, cluster_columns_louvain, "annotation")

num_pat <- length(unique(clustered_data$patient))
for (partition_size in partition_sizes){
    for (partition_index in 1:num_random_partitions){
        partition_index_str <- paste0("partition_",partition_size,"_",partition_index)
        partition_names <- paste0("Batch",1:partition_size)
        set.seed(partition_index)
        partition_assignments <- sample(x=partition_names,size=num_pat,replace=TRUE)
        names(partition_assignments) <- unique(clustered_data$patient)
        partition_vector <- sapply(clustered_data$patient, function(pat_label){partition_assignments[pat_label]})
        clustered_data[,partition_index_str] <- partition_vector
    }   
}

sample_metadata <- unique(clustered_data[,c(c("patient","sample","group","batch"),partition_name_vec)])
rownames(sample_metadata) <- sample_metadata[,"sample"]
sorted_samples <- as.character(sort(sample_metadata[,"sample"]))
sample_metadata <- sample_metadata[sorted_samples,]

options(warn=-1)
results_list <- list()

for (patition_index in 1:length(partition_name_vec)){
    partition_name <- partition_name_vec[patition_index]
    partition_name_parse <- strsplit(partition_name,"_")[[1]]
    num_partitions <- partition_name_parse[2]
    partition_replicate <- partition_name_parse[3]
    
    results_list[[partition_name]] <- list() 
    for (clustering_index in 1:length(cluster_columns)){
        clustering_name <- cluster_columns[clustering_index]
        
        cluster_assignments <- as.factor(clustered_data[,clustering_name])
        num_clusters <- length(unique(cluster_assignments))
        
        method_value <- "GloScopeProp"
        if (clustering_name == "annotation"){
            algorithm_value <- "Annotation"
            k_value <- NA
            resolution_value <- NA
            seed_value <- 1
        } else{
            col_name_parse <- strsplit(clustering_name,"_")[[1]]
            algorithm_value <- col_name_parse[1]
            k_value <- sub('.', '', col_name_parse[2])
            resolution_value <- sub('.', '', col_name_parse[3])
            seed_value <- sub('.', '', col_name_parse[4])
        }
        
        # Standardize the GS matrix and take a root transform
        gloscope_proportion_mat <- gloscope_proportion(clustered_data$sample,
            cluster_assignments, ep = 0.5, dist_mat = "KL")
        gloscope_proportion_mat <- gloscope_proportion_mat[sorted_samples,sorted_samples] # IMPORTANT!
        gloscope_proportion_mat <- apply(gloscope_proportion_mat,c(1,2),function(x){max(x,1e-16)})
        diag(gloscope_proportion_mat) <- 0
            
        # Get the separation statistics
        stats_vec <- get_statistics(gloscope_proportion_mat,
            as.integer(as.factor(sample_metadata[,partition_name])))
    
        results_vector <- c(method = method_value, algorithm = algorithm_value,
            k = k_value, resolution = resolution_value, replicate = seed_value,
            R = stats_vec$R, O2 = stats_vec$O2, silhouette = stats_vec$silhouette,
            num_clusters = num_clusters, separation_var = partition_name, gs_transform = "none",
            num_partitions = num_partitions, partition_replicate = partition_replicate)

        results_list[[partition_name]][[clustering_name]] <- results_vector
    }
}

options(warn=0)

# Parse data

results_list_byPartiton <- lapply(results_list,function(sub_list){
    sub_list_vectors <- lapply(1:length(sub_list), function(i){
        cluster_name <- names(sub_list)[i]
        result_vec <- sub_list[[i]]
        return(result_vec)
    })
    sub_list_data <- do.call(rbind, sub_list_vectors)
    sub_list_data <- as.data.frame(sub_list_data)
    return(sub_list_data)
})

results_list_final <- lapply(1:length(results_list_byPartiton),function(i){
    partition_name <- names(results_list_byPartiton)[i]
    partition_data <- results_list_byPartiton[[i]]
    return(partition_data)
})

results_df <- do.call(rbind, results_list_final)

saveRDS(results_df, here::here("results","Revision","cluspropFinal","pelkaCRC_partition_results.RDS"))