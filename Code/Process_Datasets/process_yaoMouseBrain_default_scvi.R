library(here)
library(SingleCellExperiment)
source(here::here("Code","Process_Datasets","run_scvi.R"))

# run scvi for the default Allen Mouse brain cell processing
# reduction keys must be set manually when updating the Seurat object with scVI embeddings
# see `run_scvi.R` for details about how to name reduction keys
sce_object <- readRDS(here::here("data","Processed_Datasets","yaoMouseBrain","yaoMouseBrain_default","yaoMouseBrain_default.Rds")) 

default_scvi <- run_scvi(sce_object,batch_var=NULL,is_sce=TRUE)

SingleCellExperiment::reducedDim(sce_object, "scvi", withDimnames=FALSE) <- default_scvi

saveRDS(sce_object,here::here("data","Processed_Datasets","yaoMouseBrain","yaoMouseBrain_default","yaoMouseBrain_default.Rds"))
