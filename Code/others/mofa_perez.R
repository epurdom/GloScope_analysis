library(dplyr)
library(stringr)
myseurat <- readRDS("../../data/Processed_Datasets/perezLupus/metadata/perezLupus_default_metadata.Rds")
load("../../results/BatchStudy/Perez2022_pseudobulk_manual.Rda")
meta <- unique(myseurat[, c("group", "sample")])
for(i in names(pb_cluster)){
  colnames(pb_cluster[[i]]) = paste0(colnames(pb_cluster[[i]]), "_", i)
}
pb_count <- do.call(cbind, pb_cluster)
metadata <- data.frame(id = colnames(pb_count), cell_counts = colSums(pb_count))
metadata$donor_id <- unlist(lapply(str_split(metadata$id, "_"), function(x) paste0(x[1], ".", x[2])))
metadata$cell_type <- unlist(lapply(str_split(metadata$id, "_"), function(x) x[3]))

library(MOFAcellulaR)
multiview_dat <- MOFAcellulaR::create_init_exp(counts = pb_count,  
                                               coldata = metadata)
multiview_dat <- MOFAcellulaR::filt_profiles(pb_dat = multiview_dat,
                                             cts = unique(metadata$cell_type),
                                             ncells = 1, 
                                             counts_col = "cell_counts", 
                                             ct_col = "cell_type")
multiview_dat <- MOFAcellulaR::filt_views_bysamples(pb_dat_list = multiview_dat,
                                                    nsamples = 0)

multiview_dat <- MOFAcellulaR::filt_gex_byexpr(pb_dat_list = multiview_dat,
                                               min.count = 5, # Modify!!
                                               min.prop = 0.25) # Modify!!

multiview_dat <- filt_views_bygenes(pb_dat_list = multiview_dat,
                                    ngenes = 15)
multiview_dat <- filt_samples_bycov(pb_dat_list = multiview_dat,
                                    prop_coverage = 0.9)
multiview_dat <- MOFAcellulaR::tmm_trns(pb_dat_list = multiview_dat,
                                        scale_factor = 1000000)
multiview_dat <- pb_dat2MOFA(pb_dat_list = multiview_dat, 
                             sample_column = "donor_id")
MOFAobject <- MOFA2::create_mofa(multiview_dat)

data_opts <- MOFA2::get_default_data_options(MOFAobject)
train_opts <- MOFA2::get_default_training_options(MOFAobject)
model_opts <- MOFA2::get_default_model_options(MOFAobject)

# This avoids the regularization of multicellular programs per cell type.
# This avoids less sparse gene weights
model_opts$spikeslab_weights <- FALSE 

# Define the number of factors needed
model_opts$num_factors <- 5

# Prepare MOFA model:
MOFAobject <- MOFA2::prepare_mofa(object = MOFAobject,
                                  data_options = data_opts,
                                  model_options = model_opts,
                                  training_options = train_opts)
model <- MOFA2::run_mofa(MOFAobject, save_data  = FALSE, use_basilisk  = TRUE)


save(model, file = "../../results/BatchStudy/mofa_test_perez.Rda")
