library(here)
library(BiocParallel)                                                        
source(here::here("Code","Analyze_Datasets","get_distances.R"))

parallel_params <- BiocParallel::SerialParam(RNGseed=2)  
program_vec <- rep("liuRash31_default",6)
dist_params <- list(
	list(dim_reduction="PC",sample_id="sample",ndim=10,dist="KL",
		dens="GMM",r=10000,parallel=parallel_params,
		varapp=FALSE,epapp=FALSE),
	list(dim_reduction="PC",sample_id="sample",ndim=10,dist="KL",
	     dens="KNN",k=10,parallel=parallel_params),
	list(dim_reduction="harmony",sample_id="sample",ndim=10,dist="KL",
		dens="GMM",r=10000,parallel=parallel_params,
		varapp=FALSE,epapp=FALSE),
	list(dim_reduction="harmony",sample_id="sample",ndim=10,dist="KL",
		dens="KNN",k=10,parallel=parallel_params),
	list(dim_reduction="scvi",sample_id="sample",ndim=10,dist="KL",       
                dens="GMM",r=10000,parallel=parallel_params,
                varapp=FALSE,epapp=FALSE),                                   
        list(dim_reduction="scvi",sample_id="sample",ndim=10,dist="KL",
             dens="KNN",k=10,parallel=parallel_params)
)

get_distances(program_vec,dist_params)
