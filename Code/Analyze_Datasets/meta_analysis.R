library(dplyr)
library(ggplot2)
# This file contains functions which enable the comparison of `GloS` run on
# different data sets and processing pipelines

build_grand_tibble <- function(){
	# this function recurses through the `results` folder to combine all 
	# saved pairwise distance long tibbles into a single long tibble

	# get all files in the Processed_Datasets folder that end with `dist_tibble.R`
	tibble_files <- list.files(path=here::here("results","Processed_Datasets"),recursive=TRUE,full.names=TRUE,pattern="(?i)*tibble.Rds$")
	tibble_files <- tibble_files[!grepl("/Meta/", tibble_files) ] # Ignore Metadata tibble
       	#grand_tibble <- do.call(rbind,lapply(tibble_files,function(x){ readRDS(here::here(x)) }))	
	grand_tibble <- do.call(rbind,lapply(seq_along(tibble_files),function(x){ 
				tib <- readRDS(here::here(tibble_files[[x]]));
				tib$index = factor(x);
				tib$batch_corr <- ifelse(grepl("pat",tib$dim_reduction,fixed=TRUE),"patient",ifelse(grepl("batch",tib$dim_reduction,fixed=TRUE),"batch",NA));
				return(tib) }))

	# save in the `Meta` results folder with a time stamp
	saveRDS(grand_tibble, here::here("results","Processed_Datasets","Meta",
		paste("grand_tibble_",format(Sys.time(), '%Y-%m-%d-%H-%M-%S'),".Rds",sep="")))
	return(grand_tibble)
}

build_metadata_tibble <- function(){
	metadata_files <- list.files(path=here::here("data","Processed_Datasets"),recursive=TRUE,full.names=TRUE,pattern="(?i)*default_metadata.Rds$")
	metadata_tibble <- do.call(rbind,lapply(seq_along(metadata_files),function(x){
				meta_df <- readRDS(here::here(metadata_files[[x]]));
				meta_tib <- as_tibble(meta_df);
				meta_tib <- select(meta_tib,c(dataset,patient,group,sample,cell_type,batch,num_features,num_reads))
				return(meta_tib) }))
	saveRDS(metadata_tibble,here::here("results","Processed_Datasets","Meta","metadata_tibble.Rds"))
	return(metadata_tibble)
}


generate_bar_plot <- function(dist_tibble, metric_name, experiment_vec, truncate = FALSE, filename=NULL, format = "png"){
	# A dodged box and whisker chart of distance estimates within an experiment
	colours <- c("#999999")
	plot_tibble <- dist_tibble %>%
		mutate(Method = interaction(dim_reduction, dens, sep="_")) %>%
		mutate(Experiment = interaction(dataset, program, sep="_")) %>%
		filter(((!is.na(k)) | (n > 1000))) %>% 
		filter(Experiment %in% experiment_vec)
	if (truncate){
		plot_tibble <- plot_tibble %>%
			group_by(Method, dataset, index) %>% 
			filter(distance < quantile(distance,probs=c(0.98))) %>%
			ungroup()
	}
	dodged_comparison <- plot_tibble %>%       
		filter(metric == !!metric_name) %>%		
       		ggplot(aes(x=Method,y=sqrt(distance))) +
       		geom_boxplot(aes(fill=index),position = "dodge2",width=0.1, color="black", alpha=0.2) +	
       		labs(title="Distance Distributions Across Experiments",x="Dataset",y=paste("Root ",metric_name," Distance")) +
		facet_wrap(~dataset,scales="free_x") + 
		scale_x_discrete(guide = guide_axis(angle = 90))
	if(is.null(filename)){ filename <- paste0("dataset_comparison_box_",format(Sys.time(), '%Y-%m-%d-%H-%M-%S')) }
	ggsave(here::here("results","Processed_Datasets","Meta",paste0(filename,".",format)),height=20,width=20,device=format)
	return()
}

generate_dodged_comparison <- function(distance_results_tibble, filename=NULL, format = "png"){
	dodged_comparison <- distance_results_tibble %>%
		filter(dim_reduction %in% c("PC","scvi")) %>% #,"scvibatch")) %>%
		mutate(Method = interaction(dim_reduction, dens, sep="_")) %>%
		mutate(Experiment = interaction(dataset, program, sep="_")) %>%
		filter(((!is.na(k)) | (n > 1000))) %>%
		filter(metric == "KL")
	# filter ledergor outliers
	dodged_comparison <- dodged_comparison %>%
		group_by(Method,dataset,index) %>%
		filter((Method != "PC_GMM" | dataset != "ledergor") | (distance < quantile(distance,probs=c(0.98)))) %>%
		ungroup()
	output_plot <- ggplot(dodged_comparison) +
		geom_boxplot(aes(x=dataset,y=sqrt(distance), colour = Method)) +
		labs(title="Distance Distributions Across Experiments",x="Dataset",y="Root KL Divergence")
	if(is.null(filename)){ filename <- paste0("dataset_comparison_box_",format(Sys.time(), '%Y-%m-%d-%H-%M-%S')) }
	ggsave(here::here("results","Processed_Datasets","Meta",paste0(filename,".",format)),height=20,width=20,device=format)
	return()
		
}
