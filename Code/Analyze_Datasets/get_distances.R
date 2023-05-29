library(here)
library(tibble)
library(tidyr)
library(dplyr)
library(devtools)
library(stringr)
devtools::load_all(here::here("..","popPackage"))

calc_dens = function(df_list, dens = "GMM", k = 50, num_components = c(1:9),
                     BPPARAM = BiocParallel::bpparam()){
  if(dens == "GMM"){
    mod_list <- list()
    for (i in 1:length(df_list)){
	    print(names(df_list)[i])
	    mod_list[[names(df_list)[i]]] <- densityMclust(df_list[[i]], G = NULL, verbose = TRUE, plot = FALSE)
    }
    #mod_list <- BiocParallel::bplapply(df_list, function(z) densityMclust(z, G = num_components, verbose = TRUE, plot = FALSE),
    #                          BPPARAM=BPPARAM)
  }else if(dens == "KNN"){
    mod_list <- df_list
  }
  return(mod_list)
}

#' Compute Distances from a Processed Seurat Object 
#' 
#' @param program_vec (vector of strings) The names of processing programs to consider e.g. "arnon_default"
#' @param dist_list (list of lists) The list of distance parameters associated with the program at the same index
#' @examples 
#' 
get_distances <- function(program_vec,dist_list,sce=FALSE){
	# TODO: Make inputs more safe
	# verify some data is submitted
	stopifnot(length(program_vec)>0)
	# check that there is a bijection from data to distance specification
	stopifnot(length(program_vec)==length(dist_list))

	dist_tibble_list <- list() # to be populated with long pairwise distance tibbles
	
	for (index in 1:length(program_vec)){
		# parse out the name of the data and program
		name_str <- program_vec[index]
		split_name_str <- strsplit(name_str,"_")[[1]]
		dataset_name <- split_name_str[1]
		program_name <- split_name_str[2]
		save_name <- paste0(dataset_name,"_",ifelse(endsWith(program_name,".scvi"),gsub('.{5}$', '', program_name),program_name))

		distance_params <- dist_list[[index]] # load distance parameters
		# the `gloscope` functions expects a data.frame input extracted from a Seurat object
		reduction_df <- seurat_to_df(dataset_name,program_name,save_name,distance_params$dim_reduction,sce=sce)
		embedding_matrix <- reduction_df[,stringr::str_detect(colnames(reduction_df),distance_params$dim_reduction)]
		embedding_matrix <- embedding_matrix[,1:distance_params$ndim]
		cell_sample_ids <- reduction_df$sample
		# compute the distance matrix for a given experiment
		set.seed(2)
		start_time <- proc.time()[3]
		if (distance_params$dens == "KNN"){
			distance_matrix <- gloscope(embedding_matrix, cell_sample_ids,
				dens = "KNN", dist_mat = distance_params$dist,
				k = distance_params$k, BPPARAM=distance_params$parallel)
		} else if (distance_params$dens == "GMM"){
			distance_matrix <- gloscope(embedding_matrix, cell_sample_ids,
				dens = "GMM", dist_mat = distance_params$dist,
				r = distance_params$r, num_components = distance_params$num_components,
				BPPARAM=distance_params$parallel)
		} else {
			stop("Invalid density method specified")
		}
		end_time <- proc.time()[3]
		delta_time <- end_time - start_time
		print(paste0("Runtime for ", distance_params$dens," and ",distance_params$dim_reduction," was ", delta_time))
		# Save distance matrix and parameters
		save_time <- format(Sys.time(), '%Y-%m-%d-%H-%M-%S')
		if (!dir.exists(here::here("results","Processed_Datasets",dataset_name))){ dir.create((here::here("results","Processed_Datasets",dataset_name))) }
		if (!dir.exists(here::here("results","Processed_Datasets",dataset_name,save_name))){ dir.create((here::here("results","Processed_Datasets",dataset_name,save_name))) }
		saveRDS(distance_matrix, here::here("results","Processed_Datasets",dataset_name,save_name,
			paste(name_str,save_time,"_dist_mat.Rds",sep="")))
		saveRDS(distance_params, here::here("results","Processed_Datasets",dataset_name,save_name,
			paste(name_str,save_time,"_dist_params.Rds",sep="")))

		# Convert to pairwise long tibble with distance parameters
		dist_tibble <- dist_to_tibble(distance_matrix,dataset_name,program_name,distance_params)
		dist_tibble_list[[index]] <- dist_tibble
		saveRDS(dist_tibble, here::here("results","Processed_Datasets",dataset_name,save_name,
			paste(name_str,save_time,"_dist_tibble.Rds",sep="")))
	}

	return(dist_tibble_list)
}

#' Extract dimensionality reduced coordinates as data.frame from Seurat object 
#' 
#' @param dataset_name (string) Name of the dataset
#' @param program_name (string) Name of the processing program
#' @param dim_reduction (string) Name of dimensionality reduction to extract
#' @return 
#' @examples 
#' 
seurat_to_df <- function(dataset_name,program_name,save_name,dim_reduction, sce = FALSE){
	name_str <- paste0(dataset_name,"_",program_name)
	seurat_object <- readRDS(here::here("data","Processed_Datasets",dataset_name,save_name,paste0(name_str,".Rds")))

	if (sce == FALSE){
		if (dim_reduction=="PC" | dim_reduction=="pca"){
			reduction_df <- cbind(seurat_object@meta.data, seurat_object@reductions$pca@cell.embeddings)
		} else if (dim_reduction=="harmony") {
			reduction_df <- cbind(seurat_object@meta.data, seurat_object@reductions$harmony@cell.embeddings)
		} else if (dim_reduction=="harmonybatch") {
			reduction_df <- cbind(seurat_object@meta.data, seurat_object@reductions$harmony.batch@cell.embeddings)
		} else if (dim_reduction=="harmonypat") {
			reduction_df <- cbind(seurat_object@meta.data, seurat_object@reductions$harmony.pat@cell.embeddings)
		} else if (dim_reduction=="scvi") {
			reduction_df <- cbind(seurat_object@meta.data, seurat_object@reductions$scvi@cell.embeddings)
		} else if (dim_reduction=="scvibatch") { 
			reduction_df <- cbind(seurat_object@meta.data, seurat_object@reductions$scvi.batch@cell.embeddings)
		} else if (dim_reduction=="scvipat") {
			reduction_df <- cbind(seurat_object@meta.data, seurat_object@reductions$scvi.pat@cell.embeddings)
		} 
	} else if (sce == TRUE){
		if (dim_reduction=="PC" | dim_reduction=="pca"){
			pca_matrix <- seurat_object@int_colData@listData$reducedDims$PCA
			pca_df <- as.data.frame(pca_matrix)
			meta_matrix <- seurat_object@colData
			meta_df <- as.data.frame(meta_matrix)
			reduction_df <- cbind(meta_df,pca_df)
			# TODO: Is this the right place to correct column name formatting?
			colnames(reduction_df) <- gsub("PC","PC_",colnames(reduction_df),fixed=TRUE)
		} else if (dim_reduction=="HARMONY" | dim_reduction=="harmony"){
			harmony_matrix <- seurat_object@int_colData@listData$reducedDims$HARMONY
			harmony_df <- as.data.frame(harmony_matrix)
			meta_matrix <- seurat_object@colData
			meta_df <- as.data.frame(meta_matrix)
			reduction_df <- cbind(meta_df,harmony_df)
		}	
	}
	return(reduction_df)
}

dist_to_tibble <- function(distance_matrix,dataset_name,program_name,distance_params){
	# TODO: Include distance params as covariates
	distance_tibble <- tibble::as_tibble(distance_matrix)
	# fill the diagonal and lower half of the distance matrix with NA to avoid duplicate data
	distance_tibble[lower.tri(distance_tibble,diag=T)] <- NA
	distance_tibble <- distance_tibble %>% dplyr::mutate(first_unit=as.factor(colnames(distance_tibble)))
	distance_tibble_long <- tidyr::pivot_longer(distance_tibble,cols=-c(first_unit),names_to="second_unit",names_transform=list(second_unit=as.factor),values_to="distance")
	distance_tibble_long <- tidyr::drop_na(distance_tibble_long,"distance")
	distance_tibble_long <- distance_tibble_long %>%
		dplyr::mutate(dataset = dataset_name, program = program_name,
			      dim_reduction=as.factor(distance_params$dim_reduction),
			      ndim = distance_params$ndim,
			      metric = as.factor(distance_params$dist),
			      dens = as.factor(distance_params$dens),
			      n = distance_params$n, k=distance_params$k,
			      num_components = as.factor(paste(distance_params$num_components,collapse="_")),
			      varapp = distance_params$varapp, epapp = distance_params$epapp)
	# add missing columns with NA
	all_columns <- c("dataset","program","first_unit","second_unit","distance",
			 "metric","dim_reduction","ndim","dens","n","k","num_components","varapp","epapp")
	missing_columns <- setdiff(all_columns,colnames(distance_tibble_long))
	distance_tibble_long[,missing_columns] <- NA
	#distance_tibble_long <- distance_tibble_long %>% mutate(ifelse(num_components=="",NA,num_components)) # catch edge case
	#browser()
	#distance_tibble_long <- distance_tibble_long %>% mutate_if(is.factor, na_if, y = as.factor(""))
	# order columns
	distance_tibble_long <- distance_tibble_long %>% select(all_columns)

	return(distance_tibble_long)
}
