#This script generates visualizations of the correlation between GMM and kNN as a function of cells per samples (sample depth).

library(dplyr)
library(ggplot2)
library(ggpubr)
library(here)
library(RColorBrewer)

grand_tibble <- readRDS(here::here("results","fig","publication_distance_tibble.Rds"))

build_correlation_plot <- function(grand_tibble,metric_name="KL"){
        plot_tibble <- grand_tibble %>%                                      
		filter(!(dataset %in% c("wuBreastCancer"))) %>%
                mutate(distance = ifelse(is.na(sqrt(distance)),0,sqrt(distance))) %>% # take the root distance
                group_by(dataset,dim_reduction,dens,k,n) %>% # we only keep one trial per cross of these variables, the one with the highest file index
                mutate(first_index=dplyr::first(index)) %>%                  
                ungroup() %>%                                                
                filter(index == first_index) %>%                             
                mutate(Method = interaction(dim_reduction, dens, sep="_")) %>% # our groups
                mutate(Experiment = interaction(dataset, program, sep="_"))

	pc_tibble <- plot_tibble %>% filter(dim_reduction=="PC")
	pc_tibble_baseline <- pc_tibble %>% filter(dens=="GMM") %>% mutate(baseline_distance=distance) %>% select(baseline_distance)	
	pc_tibble_alt <- pc_tibble %>% filter(dens=="KNN")
	pc_tibble <- cbind(pc_tibble_alt,pc_tibble_baseline)

	scvi_tibble <- plot_tibble %>% filter(dim_reduction=="scvi")
	scvi_tibble_baseline <- scvi_tibble %>% filter(dens=="GMM") %>% mutate(baseline_distance=distance) %>% select(baseline_distance)
	scvi_tibble_alt <- scvi_tibble %>% filter(dens=="KNN")
	scvi_tibble <- cbind(scvi_tibble_alt,scvi_tibble_baseline)

	distance_tibble <- rbind(pc_tibble,scvi_tibble) %>%
		mutate(`Dimensionality Reduction` = dim_reduction)
	levels(distance_tibble$dataset) <- c("Melanoma Tumors [Arnon et al. (2018)]",
			"Rash Epidermides [Cheng et al. (2018)]", 
			"Multiple Myeloma Blood Plasma [Ledergor et al. (2018)]",
			"Rash Immune Cells [Liu et al. (2022)]",
			"COVID Lung Tissue [Melms et al. (2021)]",
			"Colorectal Tumors [Pelka et al. (2021)]",
			"Lupus Immune Cells (Perez et al. (2022)]",
			"COVID Immune Cells [Stephenson et al. (2021)]",
			"Mouse Brain Tissue [Yao et al. (2021)]")
	levels(distance_tibble$dens) <- c("GMM","kNN")
        correlation_tibble <- distance_tibble %>%                                
                group_by(Experiment, Method, dataset,dim_reduction) %>%                             
                arrange(Experiment, Method, dataset, dim_reduction) %>%                              
                summarise(correlation = cor(distance,baseline_distance))  
	
	pca_tibble <- distance_tibble %>%
		filter(dim_reduction=="PC")
	pca_plot <- ggplot2::ggplot(pca_tibble, aes(x=baseline_distance,y=distance)) +
		geom_point(color="#009E73", alpha = 0.20) +
		geom_abline(slope = 1, intercept = 0, linetype = "longdash") +
		facet_wrap(vars(dataset), ncol=2) + #, scales="free_x") +
		theme_bw() +
		geom_smooth(method="loess",se=FALSE, linetype = "dotdash", fill = "#CC79A7") +
		ggpubr::stat_cor(aes(label = ggplot2::after_stat(r.label)), color = "black", geom = "label", cor.coef.name = "R", size=4)+
		labs(x = paste0("GMM Root Symmetric KL Divergence"),
		     y = paste0("kNN Root Symmetric KL Divergence")) +
		theme(axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 12),
		      axis.title.y = element_text(size = 24), axis.text.y = element_text(size = 12),
		      strip.text = element_text(size = 16))
                                             
        ggsave(here::here("results","fig",paste0("realworld_summary_pca.pdf")),height=20,width=20,device="pdf")
                                       

	scvi_tibble <- distance_tibble %>%
		filter(dim_reduction=="scvi")
	scvi_plot <- ggplot2::ggplot(scvi_tibble,  aes(x=baseline_distance,y=distance)) +
		geom_point(color="#009E73", alpha = 0.20) +
		geom_abline(slope = 1, intercept = 0, linetype = "longdash") +
		facet_wrap(vars(dataset), ncol=2) +
		theme_bw() +
		geom_smooth(method="loess",se=FALSE, linetype = "dotdash", fill = "#CC79A7") +
		ggpubr::stat_cor(aes(label = ggplot2::after_stat(r.label)), color = "black", geom = "label", cor.coef.name = "R", size=4)+
		labs(x = paste0("GMM Root Symmetric KL Divergence"),
		     y = paste0("kNN Root Symmetric KL Divergence")) +
		theme(axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 12),
		      axis.title.y = element_text(size = 24), axis.text.y = element_text(size = 12),
		      strip.text = element_text(size = 16))

	ggsave(here::here("results","fig",paste0("realworld_summary_scvi.pdf")),height=20,width=20,device="pdf")
        return()                                                  
}              

build_correlation_plot(grand_tibble,metric_name="KL")
