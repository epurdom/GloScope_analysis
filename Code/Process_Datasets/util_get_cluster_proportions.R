library(tidyverse)

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

get_patient_sample_cross <- function(metadata_df){
	metadata_tibble <- tibble(metadata_df)
	count_tibble <- metadata_tibble %>%
		group_by(patient,sample) %>%
		summarise(n = n())
	single_sample <- (sum(count_tibble$n) == dim(metadata_tibble)[1])
	return(single_sample)
}
