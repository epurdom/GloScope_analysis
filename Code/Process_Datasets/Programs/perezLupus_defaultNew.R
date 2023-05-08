library(here)
library(harmony)                                                             
library(SingleCellExperiment)                                                
library(scuttle) 
library(scater)
library(zellkonverter)
library(HDF5Array)

#' Format single cell data as a Seurat object
#' 
#' This function transforms raw data found in the `pop_data/data/Processed_Datasets/{dataset_name}/raw_files` folder into a `Seurat` object.
#' The `Seurat` object is minimal containing only the raw counts and original metadata, before any further processing.
#' The output Seurat object can also saved to `pop_data/data/Processed_Datasets/{dataset_name}` in a number of formats
#'
#' @param dataset_name (String) Name of folder with a dataset's raw files
#' @param assay_name (String)[optional] Name of the assay type. Default: "RNA"
#' @param save_formats (Vector of strings)[optional] Vector of other formats to save the seurat object as. Currently supported: "seurat", "list", "sce" (SingleCellExperiment), "h5Seurat". Default: c("seurat")
#' @return seurat_object (Seurat object) Processed Seurat object OR seurat_h5_filepath (String) A filepath for the connection to an h5Seurat object
#' @examples 
#' test_object <- build_seurat_object("arnon")
#' print(test_object[[]]) # get Seurat metadata

build_seurat_object <- function(dataset_name,program_name,assay_name="RNA",save_formats=c("seurat")){

	name_str <- paste0(dataset_name,"_",program_name)
        # Run checks on where the data is located
        if (!(assertthat::is.string(dataset_name) && dataset_name!= "")){ stop("Invalid dataset name") } # check the the dataset name is a valid string
        if (!dir.exists(here::here('data','Processed_Datasets',dataset_name,'raw_files'))){ stop(paste0("Folder with original data not found at:",here::here('data','Processed_Datasets','raw_files',dataset_name))) }

        #####
        # This block creates the `count_data` matrix-type object which holds the count data for the Seurat object
        # IMPORTANT: For proper interface with Seurat, this should have features as rows and cells as columns
        # Seurat has a special function for loading `.mtx` data. See: `Seurat::ReadMtx()`
        # Seurat has a special function for loading `.hdf5` data. See: `SeuratDisk::ReadH5()` which is invoked by `as.data.frame(hdf5_object)`
	# IMPORANT: Use `here::here()` for file paths relative to `pop_data`
        ####
    
	system(paste0("python3.10 ",here::here("Code","Process_Datasets","Programs","perezLupus_default.py"))) # convert raw files to readable .h5ad file

	raw_file_path <- here::here("data","Processed_Datasets","perezLupus","raw_files","scRNA_raw_cleaned.h5ad")
	log_object <- HDF5Array::H5ADMatrix(raw_file_path, "data")
	count_object <- HDF5Array::H5ADMatrix(raw_file_path, "counts")
	sce_object <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = count_object, log = log_object))

        #####
        # This block creats the associated `metadata_df` data.frame. IMPORANT: The data.frame class is required!
        # In contrast to `count_data`, cells are rows and metadata features are columns
        # The `metadata_df` object should NOT be further processed. This is done at a later stage in the pipeline
        # If no metadata is available, leave the variable assignment NULL       
        #####

	# the original object loads zero counts but has complete metadata
	original_object <- zellkonverter::readH5AD(here::here("data","Processed_Datasets","perezLupus","raw_files","local.h5ad"))
	metadata_object <- colData(original_object)
	rownames(metadata_object) = NULL                                                     
	#metadata_object <- readRDS(here::here("data","Processed_Datasets",dataset_name,"raw_files","metadata.Rds"))

	sce_object@colData <- metadata_object
	metadata_df <- as.data.frame(metadata_object)
	saveRDS(metadata_df, here::here("data","Processed_Datasets","perezLupus","raw_files","metadata_df.Rds"))                     
        if(!(is.null(metadata_df) || is.data.frame(metadata_df))){ stop("Invalid Format: Metadata is not either a data.frame or NULL") }

        #####################
        # DO NOT EDIT BELOW #
        #####################

	ifelse(!dir.exists(here::here("data","Processed_Datasets",dataset_name,name_str)), dir.create(here::here("data","Processed_Datasets",dataset_name,name_str)), FALSE)

        if ("sce" %in% save_formats){
                seurat_sce_filepath <- here::here("data","Processed_Datasets",dataset_name,name_str,"initial_sce.Rds")
                saveRDS(sce_object,file=seurat_sce_filepath)
        }   

        return(sce_object)
}


#' Custom metadata manipulation 
#' 
#' This function performs any custom manipulations of an experiment's metadata.
#' For example this could include creating new covariates from original ones or changing data types
#' In the processing pipeline this function is run on the clean Seurat object and before standardization of common covariate names
#'
#' @param metadata_df (data.frame) Metadata from a Seurat object 
#' @return metadata_df (data.frame) Updated metadata
#' @examples 
#' 
metadata_manipulation_fn <- function(metadata_df){

	metadata_df$final_sample <- paste0(metadata_df$sample_uuid, metadata_df$processing_cohort)
	return(metadata_df)
}

#' Get Seurat pipeline parameters
#'
#' This function returns a list of lists for arguments to the Seurat processing pipeline
#' Each element of the list shares its name with the relavent Seurat function
#' The value is NULL if the function should not be executed; alternatively that list element can be omitted
#' An empty list implies the default parameters and the sublist elements should otherwise correspond to Seurat function arguments
#' 
#' @return seurat_params (List of lists) Arguments to Seurat functions in the standard pipeline
#' @examples 
#' seurat_params <- list( NormalizeData = list(normalization.method="RC",scale.factor=10000) ) # only normalize data
#' 
get_seurat_pipeline_params <- function(){
	seurat_params <- list(
		NormalizeData = NULL,
		FindVariableFeatures = NULL,
		ScaleData = NULL,
		RunPCA = NULL,
		RunUMAP = NULL,
		FindNeighbors = NULL,
		FindClusters = NULL
	)

	# Verify the `seurat_params` list is properly format
        pipeline_fns <- c("NormalizeData","FindVariableFeatures","ScaleData","RunPCA","RunUMAP","FindNeighbors","FindClusters")
	for(name in names(seurat_params)[!(names(seurat_params) %in% pipeline_fns)]){
		warning(paste0("Seurat step ", name, " not in the standard Seurat pipeline. The corresponding function will not be executed."))
	}

	return(seurat_params)
}

#' Return a list of custom functions to dovetail into the Seurat pipeline
#'
#' This function returns a list of functions. The names of the list correspond to the pipeline step which should be executed directly before.
#' If a pipeline step is omitted or set to NULL, no additional function will be executed
#' For pre-processing before normalization, the element name is "Pre"
#' Custom functions should be defined below this function and should map Seurat objects (and perhaps params) to an updated Seurat object
#' The strcture of each element is as follows: prior_seurat = list(fn_name="string",params=list(a=2,b="c",...))
#' 
#' @return custom_fns_list (List of functions and NULL) Functions to execute in the Seurat pipeline
#' @examples 
#' 
get_seurat_pipeline_custom_fns <- function(){ 
	custom_fns_list <- list(
		Pre = NULL,
		NormalizeData = list(fn_name=sce_normalize,params=list()),
		ScaleData = NULL,
		RunPCA = list(fn_name=sce_pca,params=list()),
		RunUMAP = NULL,
		FindNeighbors = NULL,
		FindClusters = NULL # post-processing will come here
	)
	return(custom_fns_list)
}

###########################################
# Custom functions should be placed below #
###########################################

# For now SCE functions cannot take additional arguments 

sce_normalize <- function(sce_object){ 
        expression_mask <- (Matrix::rowSums(counts(sce_object) > 0) >= 10)
        sce_object <- sce_object[expression_mask,]                        

	sce_object <- scuttle::logNormCounts(sce_object)
	return(sce_object)
}

sce_pca <- function(sce_object){
        variable_genes <- scran::modelGeneVar(sce_object)                 
        selected_genes <- scran::getTopHVGs(variable_genes, n = 2000)     
        rowSubset(sce_object) <- selected_genes                           
        sce_object <- scater::runPCA(sce_object, scale = TRUE, subset_row = selected_genes)

	sce_object <- scater::runPCA(sce_object)
	sce_object <- harmony::RunHarmony(sce_object,group.by.vars="batch")
	return(sce_object)
}

##########################################
##########################################

#' Generate custom plots
#'
#' Plots appended to `plot_list` are called at the end of the processing pipeline
#' 
#' @param seurat_object (Seurat) Processed object
#' @return plot_list (list) List of plots to be returned a displayed with a call to `print()`
#' @examples 
#' 
get_custom_plots <- function(seurat_object){
	plot_list <- list() # If desired, generate custom plots and append to plot_list 
	return(plot_list)
}


