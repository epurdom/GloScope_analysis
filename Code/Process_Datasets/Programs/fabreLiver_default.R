# This file contains functions for processing a single dataset from raw files to a final Seurat object 
# In practice, make a copy and add your code at each TODO: comment.

library(here)
library(harmony)

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
    
        count_data <-  ReadMtx(
          mtx = here::here("data","Processed_Datasets","FabreLiver","raw_files","expression","64130e92cc37d750c0a62f14","matrix.mtx.gz"), 
          features = here::here("data","Processed_Datasets","FabreLiver","raw_files","expression","64130e92cc37d750c0a62f14","features_raw.tsv.gz"),
          cells = here::here("data","Processed_Datasets","FabreLiver","raw_files","expression","64130e92cc37d750c0a62f14","barcodes_raw.tsv.gz"))
        norm_data <-  ReadMtx(
          mtx = here::here("data","Processed_Datasets","FabreLiver","raw_files","expression","641314b01533ea248484b572","matrix__1_.mtx.gz"), 
          features = here::here("data","Processed_Datasets","FabreLiver","raw_files","expression","641314b01533ea248484b572","features__1_.tsv.gz"),
          cells = here::here("data","Processed_Datasets","FabreLiver","raw_files","expression","641314b01533ea248484b572","barcodes__1_.tsv.gz"))
        count_data <- count_data[rownames(norm_data),]
        #####
        # This block creats the associated `metadata_df` data.frame. IMPORANT: The data.frame class is required!
        # In contrast to `count_data`, cells are rows and metadata features are columns
        # The `metadata_df` object should NOT be further processed. This is done at a later stage in the pipeline
        # If no metadata is available, leave the variable assignment NULL       
        #####

        metadata_df <- read.table(here::here("data","Processed_Datasets","FabreLiver","raw_files","metadata","2022_SI_human_liver_allcells_metadata_file.txt"),header = T, sep="\t", quote = "")[-1,]
        if(!(is.null(metadata_df) || is.data.frame(metadata_df))){ stop("Invalid Format: Metadata is not either a data.frame or NULL") }

        #####################
        # DO NOT EDIT BELOW #
        #####################

        # Build and save Seurat object
        seurat_object <- Seurat::CreateSeuratObject(counts=count_data,data = norm_data,
                                                    project=name_str,assay=assay_name,meta.data=metadata_df)

        # Save Seurat object in formats as requested
	# make (if needed) a sub directory for this program
	ifelse(!dir.exists(here::here("data","Processed_Datasets",dataset_name,name_str)), dir.create(here::here("data","Processed_Datasets",dataset_name,name_str)), FALSE)

        if ("seurat" %in% save_formats){
                seurat_object_filepath <- here::here("data","Processed_Datasets",dataset_name,name_str,"initial_seurat.Rds")
                saveRDS(seurat_object,file=seurat_object_filepath)
        }   
        if ("sce" %in% save_formats){
                seurat_sce <- Seurat::as.SingleCellExperiment(seurat_object) 
                seurat_sce_filepath <- here::here("data","Processed_Datasets",dataset_name,name_str,"initial_sce.Rds")
                saveRDS(seurat_sce,file=seurat_sce_filepath)
        }   
        if("list" %in% save_formats){
                seurat_tolist <- SeuratObject::S4ToList(seurat_object)
                seurat_tolist_filepath <- here::here("data","Processed_Datasets",dataset_name,name_str,"initial_seurat_list.Rds")
                saveRDS(seurat_tolist,file=seurat_tolist_filepath)
        }   
        if ("h5Seurat" %in% save_formats){
                seurat_h5_filepath <- here::here("data","Processed_Datasets",dataset_name,name_str,"initial_seurat.h5Seurat")
                SeuratDisk::SaveH5Seurat(seurat_object,seurat_h5_filepath,overwrite=T)
        }   

        # In the future it is possible to support loom and AnnData (h5ad)

        return(seurat_object)
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

	factor_cols <- c("disease__ontology_label","paper")
	metadata_df[factor_cols] <- lapply(metadata_df[factor_cols], factor)

	metadata_df["sample"] <- metadata_df["biosample_id"]
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
#		NormalizeData = list(normalization.method="LogNormalize",scale.factor=10000,verbose=F),
		FindVariableFeatures = list(selection.method = "vst",nfeatures=2000,verbose=F),
		ScaleData = list(verbose=F),
		RunPCA = list(verbose=F)
		#		RunUMAP = list(reduction = "pca",dims=1:10,verbose=F),
#		FindNeighbors = list(verbose=F),
#		FindClusters = list(verbose=F)
	)

	# Verify the `seurat_params` list is properly format
        pipeline_fns <- c("NormalizeData","FindVariableFeatures","ScaleData","RunPCA")
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
		Pre = list(fn_name=quality_control,params=list()), 
		NormalizeData = NULL,
		ScaleData = NULL,
		RunPCA = list(fn_name=run_harmony,params=list()),
		RunUMAP = NULL,
		FindNeighbors = NULL,
		FindClusters = NULL # post-processing will come here
	)
	return(custom_fns_list)
}

###########################################
# Custom functions should be placed below #
###########################################

quality_control <- function(seurat_object){                                  
  remove_sample <- names(which(table(seurat_object$sample)<500))
  seurat_object <- seurat_object[,!(seurat_object$sample %in% remove_sample)]
  return(seurat_object)   
}

run_harmony <- function(seurat_object){ 
  seurat_object <- harmony::RunHarmony(seurat_object,"sample")
        seurat_object <- harmony::RunHarmony(seurat_object,"batch", reduction.save = "harPaper")
        return(seurat_object)
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


