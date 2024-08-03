here::i_am("Code/cluspropAnalysis/clusprop_helpers.R")

library(cluster)
library(reticulate)
library(GloScope)
library(lisi)
library(sceasy)
library(tidyverse)
library(vegan)

##############
# Statistics #
##############
                                
get_bootstrap_indices <- function(sample_metadata_df, condition, random_seed, batch = NULL){
    bootstrap_sample_ids <- list()
    sample_indices <- 1:nrow(sample_metadata_df)

    if (is.null(batch)){
        for(condition_label in unique(sample_metadata_df[,condition])){
            # Generate bootstrap samples within each condition 
            sample_ids <- sample_indices[sample_metadata_df[,condition] == condition_label]
            set.seed(random_seed)
            bootstrap_indices <- sample(sample_ids, replace = TRUE)
            # Add the bootstrap sample IDs to the list
            bootstrap_sample_ids[[condition_label]] <- bootstrap_indices
        }
    } else {
        for(batch_label in unique(sample_metadata_df[,batch])){
            for(condition_label in unique(sample_metadata_df[,condition])){
                # Generate bootstrap samples within each cross of condition and batch
                sample_ids <- (sample_indices[sample_metadata_df[,batch] == batch_label &
                    sample_metadata_df[,condition] == condition_label])
                set.seed(random_seed)
                bootstrap_indices <- sample(sample_ids, replace = TRUE)
                # Add the bootstrap sample IDs to the list
                bootstrap_sample_ids[[paste0(batch_label,condition_label)]] <- bootstrap_indices
            }
        }
    }
    
    bootstrap_sample_ids <- unlist(bootstrap_sample_ids)
    return(bootstrap_sample_ids)
}

get_statistics <- function(dist_mat, group_id_vec){

    anosim_result <- vegan::anosim(dist_mat, group_id_vec)
    anosim_statistic <- anosim_result$statistic

    permanova_result <- vegan::adonis2(dist_mat ~ group_id_vec)
    permanova_statistc <- get_O2_statistic(permanova_result)

    silhouette_statistic <- get_silhouette_statistic(dist_mat, group_id_vec)
    
    

    return_list <- list(R = anosim_statistic,
        O2 = permanova_statistc,
        silhouette = silhouette_statistic)

    return(return_list)
}

get_O2_statistic <- function(permanova_result){
    sumsq <- permanova_result$SumOfSqs
    n <- permanova_result$Df
    O2_statistic <- (sumsq[1]-n[1]*(sumsq[2]/n[2]))/(sumsq[3] + sumsq[2]/n[2])
    return(O2_statistic)
}

get_silhouette_statistic <- function(dist_mat, group_covar_vec){
    silhouette_result <- cluster::silhouette(x = as.integer(group_covar_vec), dmatrix = dist_mat)
    # TODO: Is this the right data transformation?
    silhouette_statistic <- mean(silhouette_result[,"sil_width"])
    return(silhouette_statistic)
}

##############################
# Analysis and Visualization #
##############################

get_cluster_props <- function(metadata_df,unit_col,label_col){
	# unit_col could be sample or patient ID, depending on preference
	metadata_tibble <- tibble(metadata_df)
	proportion_tibble <- metadata_tibble %>%
		mutate(across(c(.data[[unit_col]],.data[[label_col]]),factor)) %>%
		group_by(.data[[unit_col]],.data[[label_col]]) %>%
		summarise(n = n()) %>%
		mutate(freq = n / sum(n)) %>%
		select(-n) %>% # need to remove `n` to pivot correctly
		pivot_wider(names_from=.data[[label_col]],values_from=freq,values_fill=0)
	proportion_df <- data.frame(proportion_tibble)

	# move `unit_col` values to row names 
	rownames(proportion_df) <- proportion_df[,1]
	proportion_df <- proportion_df[,-1, drop=FALSE]
	return(proportion_df)
}
