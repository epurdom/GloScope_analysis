source("Code/Process_Datasets/run_scvi.R")
library(here)

# run scvi for the default Pelka colorectal cancer processing
# reduction keys must be set manually when updating the Seurat object with scVI embeddings
# see `run_scvi.R` for details about how to name reduction keys
seurat_object <- readRDS(here::here("data","Processed_Datasets","pelkaCRC","pelkaCRC_default","pelkaCRC_default.Rds")) 

default_scvi <- run_scvi(seurat_object,batch_var=NULL)
seurat_object[["scvi"]] <- CreateDimReducObject(embeddings = default_scvi, key = "scvi_", assay = DefaultAssay(seurat_object))

saveRDS(seurat_object,here::here("data","Processed_Datasets","pelkaCRC","pelkaCRC_default","pelkaCRC_default.Rds"))
