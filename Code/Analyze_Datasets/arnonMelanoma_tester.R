library(here)
library(BiocParallel)
source(here::here("Code","Analyze_Datasets","get_distances.R"))

program_vec <- rep("arnonMelanoma_default",1)
parallel_params <- BiocParallel::SerialParam(RNGseed=2)
dist_params <- list(
	list(dim_reduction="PC",sample_id="sample",ndim=10,dist="KL",
		dens="GMM",r=10000,parallel=parallel_params,
		varapp=FALSE,epapp=FALSE)
)
get_distances(program_vec,dist_params)
