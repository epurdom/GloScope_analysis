library(Seurat)

# This function turns count matrices into seurat objects 
# It expects a cell x gene matrix with named rows and columns, as well as the relavent metadata
# The list `seurat_params` allows the user to use non-default parameters
# TODO: Details on the format of `seurat_params`
# NULL entires in `seurat_params` will be skilled e.g. `seurat_params$umap <- NULL` means UMAP will not be run
# IMPORANT: This may break downstream analysis e.g. UMAP uses PCA by default
# `additional_steps` is a function for further processing eg tSNE or kNN. It expeects and returns a Seurat object 

# Credit: https://github.com/satijalab/seurat/issues/5343
ClearSeuratCommands <- function(seuratObj, maxSize = 5000) {
	for (commandName in names(seuratObj@commands)) {
		val <- object.size(x = slot(seuratObj@commands[[commandName]], 'call.string'))
		if (val > maxSize) {
			print(paste0('Clearing call.string for: ', commandName, '. size: ', format(val, units = 'auto')))
			slot(seuratObj@commands[[commandName]], 'call.string') <- ''
		}
	}
	return(seuratObj)
} 

#' Run a Seurat pipeline program
#'
#' This function runs a Seurat pipeline on an input object. Standard functions include normalization, finding variable features, scaling, PCA, UMAP, and clustering.
#' Arguments to each of these steps are passed in through the `seurat_params` list. Custom functions to execute after each step are passed in through the `custom_fns` list.
#' 
#' @param seurat_object (Seurat object) Input Seurat object with desired metadata
#' @param seurat_params (List of lists) List of parameters for the Seurat functions in the pipeline
#' @param custom_fns value (List of lists) List of custom functions to execute after each pipeline step
#' @return seurat_object (Seurat object) The updated Seurat object
#' @examples 
#' 
run_seurat_pipeline <- function(seurat_object,seurat_params,custom_fns=list(),is_sce=FALSE){
	custom_function_locations <- names(custom_fns)

	# If user specifies specific pre-processing, run that here
	if (("Pre" %in% custom_function_locations) && !(is.null(custom_fns$Pre))){
		if(is.function(custom_fns$Pre$fn_name)){
            if (!is_sce){
                seurat_object <- do.call(custom_fns$Pre$fn_name,c(seurat_object,custom_fns$Pre$params))
                seurat_object <- ClearSeuratCommands(seurat_object)
            } else if (is_sce){
                # SCE is incompatible with do.call apparently. TODO: Allow for optional argument
                seurat_object <- custom_fns$Pre$fn_name(seurat_object) 
            }
		}
	}
	
	# Normalize data
	if (!is.null(seurat_params$NormalizeData)) {
		seurat_object <- do.call(Seurat::NormalizeData,c(list(object=seurat_object),seurat_params$NormalizeData))
		seurat_object <- ClearSeuratCommands(seurat_object)
	}
	if (("NormalizeData" %in% custom_function_locations) && !(is.null(custom_fns$NormalizeData))){
		if(is.function(custom_fns$NormalizeData$fn_name)){
            if (!is_sce){
                seurat_object <- do.call(custom_fns$NormalizeData$fn_name,c(seurat_object,custom_fns$NormalizeData$params))
                seurat_object <- ClearSeuratCommands(seurat_object) 
            } else if (is_sce){
                seurat_object <- custom_fns$NormalizeData$fn_name(seurat_object) 
            }
		}
	}

	# Select variable features
	if (!is.null(seurat_params$FindVariableFeatures)){
		seurat_object <- do.call(Seurat::FindVariableFeatures,c(list(object=seurat_object),seurat_params$FindVariableFeatures))
		seurat_object <- ClearSeuratCommands(seurat_object)
	}
	if (("FindVariableFeatures" %in% custom_function_locations) && !(is.null(custom_fns$FindVariableFeatures))){
		if(is.function(custom_fns$FindVariableFeatures$fn_name)){
            if (!is_sce){
                seurat_object <- do.call(custom_fns$FindVariableFeatures$fn_name,
                                         c(seurat_object,custom_fns$FindVariableFeatures$params))
                seurat_object <- ClearSeuratCommands(seurat_object)
            }  else if (is_sce){
                seurat_object <- custom_fns$FindVariableFeatures$fn_name(seurat_object) 
            }
		}
	}

	# Center and scale data
	if (!is.null(seurat_params$ScaleData)){
		# Due to the aforementioned pipeline bug, custom arguments cannot be used for scaling right now
		seurat_object <- Seurat::ScaleData(seurat_object) 
		seurat_object <- ClearSeuratCommands(seurat_object)
	}
	if (("ScaleData" %in% custom_function_locations) && !(is.null(custom_fns$ScaleData))){
		if(is.function(custom_fns$ScaleData$fn_name)){
            if (!is_sce){
                seurat_object <- do.call(custom_fns$ScaleData$fn_name,c(seurat_object,custom_fns$ScaleData$params))
                seurat_object <- ClearSeuratCommands(seurat_object)
            } else if (is_sce){
                seurat_object <- custom_fns$ScaleData$fn_name(seurat_object) 
            }
		}
	}

	# Run PCA
	if (!is.null(seurat_params$RunPCA)){
		seurat_object <- do.call(Seurat::RunPCA,c(list(object=seurat_object),seurat_params$RunPCA))
		seurat_object <- ClearSeuratCommands(seurat_object)
	} 
    if (("RunPCA" %in% custom_function_locations) && !(is.null(custom_fns$RunPCA))){
            if(is.function(custom_fns$RunPCA$fn_name)){
                if (!is_sce){
                    seurat_object <- do.call(custom_fns$RunPCA$fn_name,c(seurat_object,custom_fns$RunPCA$params))
                    seurat_object <- ClearSeuratCommands(seurat_object)
                } else if (is_sce){
                    seurat_object <- custom_fns$RunPCA$fn_name(seurat_object) 
                }
            }   
    }   
	
	# Run UMAP
	if (!is.null(seurat_params$umap)){
		seurat_object <- do.call(Seurat::RunUMAP,c(list(object=seurat_object),seurat_params$RunUMAP))
		seurat_object <- ClearSeuratCommands(seurat_object)
	}
    if (("RunUMAP" %in% custom_function_locations) && !(is.null(custom_fns$RunUMAP))){
            if(is.function(custom_fns$RunUMAP$fn_name)){
                if (!is_sce){
                    seurat_object <- do.call(custom_fns$RunUMAP$fn_name,c(seurat_object,custom_fns$RunUMAP$params))
                    seurat_object <- ClearSeuratCommands(seurat_object)
                } else if (is_sce){
                    seurat_object <- custom_fns$RunUMAP$fn_name(seurat_object) 
                }
            }
    }

	# Find Neighbors
	if (!is.null(seurat_params$FindNeighbors)){
		seurat_object <- do.call(Seurat::FindNeighbors,c(list(object=seurat_object),seurat_params$FindNeighbors))
		seurat_object <- ClearSeuratCommands(seurat_object)
	} 
	if (("FindNeighbors" %in% custom_function_locations) && !(is.null(custom_fns$FindNeighbors))){
		if(is.function(custom_fns$FindNeighbors$fn_name)){
            if (!is_sce){
                seurat_object <- do.call(custom_fns$FindNeighbors$fn_name,c(seurat_object,custom_fns$FindNeighbors$params))
                seurat_object <- ClearSeuratCommands(seurat_object)
            } else if (is_sce){
                seurat_object <- custom_fns$FindNeighbors$fn_name(seurat_object) 
            }
		}
	}

	# Run Clustering
	if (!is.null(seurat_params$FindClusters)){
		seurat_object <- do.call(Seurat::FindClusters,c(list(object=seurat_object),seurat_params$FindClusters))
		seurat_object <- ClearSeuratCommands(seurat_object)
	}
    if (("FindClusters" %in% custom_function_locations) && !(is.null(custom_fns$FindClusters))){
        if(is.function(custom_fns$FindClusters$fn_name)){
            if (!is_sce){
                seurat_object <- do.call(custom_fns$FindClusters$fn_name,c(seurat_object,custom_fns$FindClusters$params))
                seurat_object <- ClearSeuratCommands(seurat_object)
            } else if (is_sce){
                seurat_object <- custom_fns$FindNeighbors$fn_name(seurat_object) 
            }
        }  
	}

	return(seurat_object)
}
