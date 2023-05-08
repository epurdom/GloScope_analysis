library(data.table)

#' Manipulates a Seurat metadata data.frame to have consistent format and apply custom transformations
#' 
#' There are a handful of covariates which we want for every dataset in a standardized format. These are the bulk of the arguments to this function.
#' Arguments such as `dataset` should be set to the name of the column which contains contains that information and NULL if it does not exits.
#' Custom steps to manipulate the metadata_df can be specified through the `custom_steps` function
#' The datatype of the standard covariates will be changed to factor or numeric, but remaining must be converted by the user in the custom_steps function 
#'
#' @param metadata_df (data.frame) The metadata of a Seurat object
#' @param dataset (String) The name of the data set where this data comes from
#' @param patient (String) The column in metadata_df which contains patient ID
#' @param group (String) The column which contains group membership/phenotype
#' @param sample (String) The column which contains the sample ID within a patient
#' @param batch (String) The column which contains an indicator of a possible batch covariate
#' @param num_features (String) The column which contains the number of features per cell
#' @param num_reads (String) The column which contains the number of reads per cell
#' @param custom_steps (Function or NULL) Any additional post-processing
#' @return metadata_df (data.frame) Updated metadata
#' @examples 
#' 
format_metadata_df <- function(metadata_df,dataset=NULL,patient=NULL,group=NULL,sample=NULL,cell_type=NULL,batch=NULL,num_features=NULL,num_reads=NULL,custom_steps=NULL){
	# TODO: The column name arguments could probably be refactored as a list

	# if the user has specified a custom function apply it here
	# this could include e.g. extracting the patient ID from a more complex string
	if ((!is.null(custom_steps) && is.function(custom_steps))){
		metadata_df <- custom_steps(metadata_df)
	}

	all_covariates <- colnames(metadata_df) # the covariates in the input data.frame
	standard_covariates <- c("patient","group","sample","cell_type","batch","num_features","num_reads") # the standard names shared across datasets
	standard_covariates_matches <- unlist( lapply(standard_covariates, function(x) { eval(parse(text = x)) }) ) # input data.frame names for standard covariates which exist
	sc_match_indices <- unlist( lapply(standard_covariates, function(x) {!is.null(eval(parse(text = x)))} ) ) # The indices in `standard_covariates` which have a matched covariate
	standard_covariates_present <- standard_covariates[sc_match_indices] # the names of `standard_covariates` which have a matched input covariate
	sc_null_indices <- !sc_match_indices # SC indices which are unmatched
	null_covariates <- standard_covariates[sc_null_indices]
	other_covariates <- all_covariates[!(all_covariates %in% standard_covariates_matches)] # these are covariates unique to the dataset at hand

	# The `setnames` function copies in place, this next step is to prevent bugs
	copied_metadata_df <- data.frame(metadata_df)
	metadata_df <- copied_metadata_df
	if(!is.null(standard_covariates_matches)){
		data.table::setnames(metadata_df,old=standard_covariates_matches,new=standard_covariates_present) # this is in-place!
	}

	# standardize data types for standard covariates
	factor_covariates <- c("dataset","patient","group","sample","cell_type","batch")
	numeric_covaraites <- c("num_features","num_reads")

	factor_cols <- factor_covariates[factor_covariates %in% standard_covariates_present] 
	if (length(factor_cols)>0){
		metadata_df[factor_cols] <- lapply(metadata_df[factor_cols], factor)
	}
	numeric_cols <- numeric_covaraites[numeric_covaraites %in% standard_covariates_present]
	if (length(numeric_cols)>0){
		metadata_df[numeric_cols] <- lapply(metadata_df[numeric_cols], as.numeric)
	}

	# add NAs
	metadata_df[,null_covariates] <- NA
	metadata_df[,"dataset"] <- as.factor(dataset)

	# reorder columns with standard covariates in the lead and dataset specific one following
	metadata_df <- metadata_df[,c("dataset",standard_covariates,other_covariates)] # reorder with shared covariates first

	print("Metadata covariates:")
	print(head(metadata_df))
	

	return(metadata_df)
}
