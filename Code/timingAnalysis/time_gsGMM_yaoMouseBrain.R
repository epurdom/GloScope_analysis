here::i_am("Code/timingAnalysis/time_gsGMM_yaoMouseBrain.R")

library(GloScope)

sce_object <- readRDS(here::here("data","Processed_Datasets","yaoMouseBrain","yaoMouseBrain_default","yaoMouseBrain_default.Rds"))

embedding_matrix <- SingleCellExperiment::reducedDim(sce_object,"PCA")[,1:10]
cell_sample_ids <- SingleCellExperiment::colData(sce_object)$sample

dens <- "GMM"
dist_mat <- "KL"

cell_sample_ids <- as.character(cell_sample_ids)
unique_sample_ids <- unique(cell_sample_ids)
names(unique_sample_ids) <- unique_sample_ids

BPPARAM <- BiocParallel::SerialParam()

density_start <- proc.time()[3]

sample_matrix_list <- lapply(unique_sample_ids,
    function(x){embedding_matrix[(cell_sample_ids==x),,drop=FALSE]})
mod_list <- .calc_dens(sample_matrix_list, dens = dens, k = NULL, BPPARAM = BPPARAM, num_components = c(1:9))

density_end <- proc.time()[3]
density_time <- density_end - density_start

divergence_start <- proc.time()[3]
sample_pairs <- utils::combn(unique_sample_ids, 2)
# Convert patient pairs to a list for BiocParallel::bplapply
patient_pair_list <- lapply(seq_len(ncol(sample_pairs)), function(i) sample_pairs[,i])
# IMPORTANT: There are additional algorithms for divergence estimation
# implemented in this package which are not accessible from the `gloscope`
# function. The optional arguments `varapp`, `epapp`, and `ep` must be manually
# set below. See `R/.calc_dist.R` for their details.
divergence_list <- BiocParallel::bplapply(patient_pair_list,
    function(w){ .calc_dist(mod_list = mod_list, s1 = w[1], s2 = w[2],
        df_list = sample_matrix_list, dist_mat = dist_mat, dens = dens,
        r = 10000, k = NULL,
        varapp = FALSE, epapp = FALSE, ep = NA)},BPPARAM=BPPARAM)

# Convert pair-wise distances to a symmetric distance matrix
divergence_vec <- unlist(divergence_list)
divergence_matrix <- matrix(0, ncol = length(unique_sample_ids),
                            nrow = length(unique_sample_ids))
rownames(divergence_matrix) <- unique_sample_ids
colnames(divergence_matrix) <- unique_sample_ids
divergence_matrix[lower.tri(divergence_matrix, diag=FALSE)] <- divergence_vec
divergence_matrix[upper.tri(divergence_matrix)] <- t(divergence_matrix)[upper.tri(divergence_matrix)]


divergence_end <- proc.time()[3]
divergence_time <- divergence_end - divergence_start
                            
timing_vector <- c("density" = density_time, "divergence" = divergence_time)
saveRDS(timing_vector,here::here("Code","timingAnalysis","timing_results","yaoMouseBrain_gsGMM.RDS"))