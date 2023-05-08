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
    
	mtx_files <- list.files(path=here::here('data','Processed_Datasets','waickmanDengue','raw_files'),pattern="*.mtx.gz",full.names=T)
	load_list <- lapply(mtx_files,parse_mtx_files)

        count_data <- do.call("cbind", lapply(load_list,function(x){ x$count_matrix }) )

        #####
        # This block creats the associated `metadata_df` data.frame. IMPORANT: The data.frame class is required!
        # In contrast to `count_data`, cells are rows and metadata features are columns
        # The `metadata_df` object should NOT be further processed. This is done at a later stage in the pipeline
        # If no metadata is available, leave the variable assignment NULL       
        #####

        metadata_df <- do.call("rbind", lapply(load_list,function(x){ x$metadata }) ) 
        if(!(is.null(metadata_df) || is.data.frame(metadata_df))){ stop("Invalid Format: Metadata is not either a data.frame or NULL") }

        #####################
        # DO NOT EDIT BELOW #
        #####################

        # Build and save Seurat object
        seurat_object <- Seurat::CreateSeuratObject(counts=count_data,project=name_str,assay=assay_name,meta.data=metadata_df)

	set.seed(2023)
	subset_ngenes <- 3000
	subset_genes <- sample(1:dim(seurat_object)[1],subset_ngenes,replace=FALSE)
	subset_ncells <- 7000
	subset_cells <- sample(1:dim(seurat_object)[2],subset_ncells,replace=FALSE)
	seurat_object <- seurat_object[subset_genes,subset_cells]

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


# This helper function takes in the filepath to a `mtx` files
# It returns the loaded `mtx` file and a metadata data.frame parsed from the filepath
parse_mtx_files <- function(mtx_filepath){
	split_path <- strsplit(mtx_filepath,c("_|/"))[[1]] 
	gsm_id <- split_path[14]
	patient_id <- split_path[15] 
	timepoint <- split_path[16] 

	gene_filepath <- here::here("data","Processed_Datasets","waickmanDengue","raw_files",paste(gsm_id,patient_id,timepoint,"genes.tsv.gz",sep="_"))
	cell_filepath <- here::here("data","Processed_Datasets","waickmanDengue","raw_files",paste(gsm_id,patient_id,timepoint,"barcodes.tsv.gz",sep="_"))

	count_matrix <- Seurat::ReadMtx(mtx_filepath,cell_filepath,gene_filepath)
	metadata_df <- data.frame(cell=colnames(count_matrix),patient=rep(patient_id,ncol(count_matrix)),timepoint=rep(timepoint,ncol(count_matrix)),row.names="cell")

	# add a suffix to cell names to prevent duplication with other units
	colnames(count_matrix)<-paste(colnames(count_matrix),gsm_id,sep="_")
	rownames(metadata_df)<-paste(rownames(metadata_df),gsm_id,sep="_")

	return_list <- list(count_matrix=count_matrix, metadata=metadata_df)
	return(return_list)
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

	metadata_df["group"] <- ifelse(grepl("Natural",metadata_df$patient),"Natural","Experimental")
	metadata_df["sample"] <- do.call(paste,c(metadata_df[c("patient","timepoint")],sep="_")) 

	factor_cols <- c("patient","timepoint","group","sample")
	metadata_df[factor_cols] <- lapply(metadata_df[factor_cols], factor)

	metadata_df["counts"] <- metadata_df["nCount_RNA"] # redefine to avoid issues with Seurat::PercentageFeatureSet

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
		RunPCA = list(verbose=F),
		RunUMAP = list(reduction = "pca",dims=1:10,verbose=F),
		FindNeighbors = list(verbose=F),
		FindClusters = list(verbose=F)
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
		Pre = list(fn_name=seurat_preprocessing,params=list()), # demonstration fn, returns unmodified input unless uncommented
		NormalizeData = list(fn_name=run_sct,params=list()),
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

seurat_preprocessing <- function(seurat_object){ 
	# steps shared in correspondence with authors
	seurat_object$percent.mt <- Seurat::PercentageFeatureSet(object = seurat_object, pattern = "^MT-")
	#seurat_object$percent.rbc <- Seurat::PercentageFeatureSet(object = seurat_object, features = c("HBA1","HBA2","HBB"))
	seurat_object <- subset(x = seurat_object, subset = num_features> 200 & num_features < 6000 & percent.mt < 10) #& percent.rbc < 5)

	expression_mask <- Matrix::rowSums(seurat_object@assays$RNA@counts > 0)
	seurat_object <- seurat_object[expression_mask >= 10,]

	return(seurat_object)
}

run_sct <- function(seurat_object){
	seurat_object <- Seurat::SCTransform(seurat_object)	
	return(seurat_object)
}

run_harmony <- function(seurat_object){
	seurat_object <- harmony::RunHarmony(seurat_object,"patient",reduction.save="harmony.pat")
	seurat_object <- Seurat::ProjectDim(seurat_object,reduction="harmony.pat",overwrite=T,verbose=F)
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


