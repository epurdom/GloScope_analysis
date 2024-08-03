here::i_am("Code/cluspropAnalysis/score_gloscopeVarpar.R")

source(here::here("Code","cluspropAnalysis","clusprop_helpers.R"))

parameter_settings <- c("G5","G10","G15","K10","K25","K50","K100")

get_gloscope_stats <- function(dataset_name, parameter_settings, metadata_sample, separation_var){
    results_list <- list()
    for (parameter_setting in parameter_settings){
        # Parse parameters
        parameter_name = substr(parameter_setting,1,1)
        parameter_value = substr(parameter_setting,2,nchar(parameter_setting))
        dist_mat_candidate = get(paste0(dataset_name,"_",parameter_setting))
        # Compute stats
        stat_vec <- unlist(get_statistics(dist_mat_candidate, metadata_sample[,separation_var]))
        # Format metadata
        if (parameter_name == "G"){
            meta_vec <- c(dataset = dataset_name, method = "GloScopeGMM", separation_var = separation_var,
                num_components = parameter_value, num_neighbors = NA)
        } else if (parameter_name == "K"){
            meta_vec <- c(dataset = dataset_name, method = "GloScopeKNN", separation_var,
                num_components = NA, num_neighbors = parameter_value)
        } else{stop()}
        # Return results
        results_vec <- c(meta_vec, stat_vec)
        results_list[[parameter_setting]] <- results_vec
    }
    results_list <- unname(results_list)
    return(results_list)   
}

# Fabre Liver

fabreLiver_metadata_cell <- readRDS(here::here("results","Revision","cluspropFinal","fabreLiver_clustered_revised.RDS"))
fabreLiver_metadata_sample <- unique(fabreLiver_metadata_cell[,c("dataset","patient","sample","group","batch")])
rownames(fabreLiver_metadata_sample) <- fabreLiver_metadata_sample$sample
sample_order <- rownames(fabreLiver_metadata_sample)

load(here::here("results","BatchStudy","dist_mat_fabre_liver_changepar.Rda"))
fabreLiver_G5 <- dist_mat_G5[sample_order,sample_order]
fabreLiver_G10 <- dist_mat_G10[sample_order,sample_order]
fabreLiver_G15 <- dist_mat_G15[sample_order,sample_order]
fabreLiver_K10 <- dist_mat_K10[sample_order,sample_order]
fabreLiver_K25 <- dist_mat_K25[sample_order,sample_order]
fabreLiver_K50 <- dist_mat_K50[sample_order,sample_order]
fabreLiver_K100 <- dist_mat_K100[sample_order,sample_order]
rm(list=c("dist_mat_G5","dist_mat_G10","dist_mat_G15","dist_mat_K10","dist_mat_K25","dist_mat_K50","dist_mat_K100"))

results_list_fabreLiver_batch <- get_gloscope_stats("fabreLiver", parameter_settings, fabreLiver_metadata_sample, "batch")
results_list_fabreLiver_group <- get_gloscope_stats("fabreLiver", parameter_settings, fabreLiver_metadata_sample, "group")
results_fabreLiver <- do.call(rbind,c(results_list_fabreLiver_batch, results_list_fabreLiver_group))

# Fabre Lung

fabreLung_metadata_cell <- readRDS(here::here("results","Revision","cluspropFinal",
    "fabreLung_clustered_revised.RDS"))
fabreLung_metadata_sample <- unique(fabreLung_metadata_cell[,c("dataset","patient","sample","group","batch")])
rownames(fabreLung_metadata_sample) <- fabreLung_metadata_sample$sample
sample_order <- rownames(fabreLung_metadata_sample)

load(here::here("results","BatchStudy","dist_mat_fabre_lung_pcachangepar.Rda"))
fabreLung_G5 <- dist_mat_G5[sample_order,sample_order]
fabreLung_G10 <- dist_mat_G10[sample_order,sample_order]
fabreLung_G15 <- dist_mat_G15[sample_order,sample_order]
fabreLung_K10 <- dist_mat_K10[sample_order,sample_order]
fabreLung_K25 <- dist_mat_K25[sample_order,sample_order]
fabreLung_K50 <- dist_mat_K50[sample_order,sample_order]
fabreLung_K100 <- dist_mat_K100[sample_order,sample_order]
rm(list=c("dist_mat_G5","dist_mat_G10","dist_mat_G15","dist_mat_K10","dist_mat_K25","dist_mat_K50","dist_mat_K100"))

results_list_fabreLung_batch <- get_gloscope_stats("fabreLung", parameter_settings, fabreLung_metadata_sample, "batch")
results_list_fabreLung_group <- get_gloscope_stats("fabreLung", parameter_settings, fabreLung_metadata_sample, "group")
results_fabreLung <- do.call(rbind,c(results_list_fabreLung_batch, results_list_fabreLung_group))

# Perez Lupus

perezLupus_metadata_cell <- readRDS(here::here("results","Revision","cluspropFinal","perezLupus_clustered_revised.RDS"))
perezLupus_metadata_sample <- unique(perezLupus_metadata_cell[,c("dataset","patient","sample","group","batch")])
rownames(perezLupus_metadata_sample) <- perezLupus_metadata_sample$sample
sample_order <- rownames(perezLupus_metadata_sample)

load(here::here("results","BatchStudy","distmat_Perez2022_changepar.Rda"))
perezLupus_G5 <- dist_mat_G5[sample_order,sample_order]
perezLupus_G10 <- dist_mat_G10[sample_order,sample_order]
perezLupus_G15 <- dist_mat_G15[sample_order,sample_order]
perezLupus_K10 <- dist_mat_k10[sample_order,sample_order]
perezLupus_K25 <- dist_mat_k25[sample_order,sample_order]
perezLupus_K50 <- dist_mat_k50[sample_order,sample_order]
perezLupus_K100 <- dist_mat_k100[sample_order,sample_order]
rm(list=c("dist_mat_G5","dist_mat_G10","dist_mat_G15","dist_mat_k10","dist_mat_k25","dist_mat_k50","dist_mat_k100"))

results_list_perezLupus_batch <- get_gloscope_stats("perezLupus", parameter_settings, perezLupus_metadata_sample, "batch")
results_list_perezLupus_group <- get_gloscope_stats("perezLupus", parameter_settings, perezLupus_metadata_sample, "group")
results_perezLupus <- do.call(rbind,c(results_list_perezLupus_batch, results_list_perezLupus_group))

# Stephenson COVID PBMC

stephensonCOVIDPBMC_metadata_cell <- readRDS(here::here("results","Revision","cluspropFinal","stephensonCOVIDPBMC_clustered_revised.RDS"))
stephensonCOVIDPBMC_metadata_sample <- unique(stephensonCOVIDPBMC_metadata_cell[,c("dataset","patient","sample","group","batch","Status")])
rownames(stephensonCOVIDPBMC_metadata_sample) <- stephensonCOVIDPBMC_metadata_sample$sample#
# Filter 143 samples down to 126
stephensonCOVIDPBMC_metadata_sample <- stephensonCOVIDPBMC_metadata_sample[stephensonCOVIDPBMC_metadata_sample$Status %in% c("Covid","Healthy"),]
sample_order <- rownames(stephensonCOVIDPBMC_metadata_sample)

load(here::here("results","BatchStudy","distmat_COVID143_changepar.Rda"))
load(here::here("results","BatchStudy","distmat_COVID143_changepar_KNN.Rda"))
stephensonCOVIDPBMC_G5 <- dist_mat_G5[sample_order,sample_order]
stephensonCOVIDPBMC_G10 <- dist_mat_G10[sample_order,sample_order]
stephensonCOVIDPBMC_G15 <- dist_mat_G15[sample_order,sample_order]
stephensonCOVIDPBMC_K10 <- apply(dist_mat_k10,c(1,2),function(x){max(x,1e-12)})[sample_order,sample_order]
stephensonCOVIDPBMC_K25 <- apply(dist_mat_k25,c(1,2),function(x){max(x,1e-12)})[sample_order,sample_order]
stephensonCOVIDPBMC_K50 <- apply(dist_mat_k50,c(1,2),function(x){max(x,1e-12)})[sample_order,sample_order]
stephensonCOVIDPBMC_K100 <- apply(dist_mat_k100,c(1,2),function(x){max(x,1e-12)})[sample_order,sample_order]
rm(list=c("dist_mat_G5","dist_mat_G10","dist_mat_G15","dist_mat_k10","dist_mat_k25","dist_mat_k50","dist_mat_k100"))

results_list_stephensonCOVIDPBMC_batch <- get_gloscope_stats("stephensonCOVIDPBMC", parameter_settings, stephensonCOVIDPBMC_metadata_sample, "batch")
results_list_stephensonCOVIDPBMC_group <- get_gloscope_stats("stephensonCOVIDPBMC", parameter_settings, stephensonCOVIDPBMC_metadata_sample, "group")
results_stephensonCOVIDPBMC <- do.call(rbind,c(results_list_stephensonCOVIDPBMC_batch, results_list_stephensonCOVIDPBMC_group))

# For consistency with de novo clustering experiments
num_results <- nrow(results_stephensonCOVIDPBMC)
dataset_vec <- c(rep("stephensonCOVIDPBMCfilter", num_results), rep("stephensonCOVIDPBMCsubset", num_results))
results_stephensonCOVIDPBMC <- rbind(results_stephensonCOVIDPBMC,results_stephensonCOVIDPBMC)
results_stephensonCOVIDPBMC[,"dataset"] <- dataset_vec

gloscope_results <- rbind(results_fabreLiver, results_fabreLung, results_perezLupus, results_stephensonCOVIDPBMC)
saveRDS(gloscope_results, here::here("results","Revision","cluspropFinal","gloscope_varparam_results.RDS"))