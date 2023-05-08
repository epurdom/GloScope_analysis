source("Code/Process_Datasets/run_scvi.R")
library(here)

# run scvi for the default ledergor processing
# reduction keys must be set manually when updating the Seurat object
# see `run_scvi.R` for details about how to name reduction keys
seurat_object <- readRDS(here::here("data","Processed_Datasets","ledergorMyeloma","ledergorMyeloma_default","ledergorMyeloma_default.Rds") )

default_scvi <- run_scvi(seurat_object,batch_var=NULL)
seurat_object[["scvi"]] <- CreateDimReducObject(embeddings = default_scvi, key = "scvi_", assay = DefaultAssay(seurat_object))

saveRDS(seurat_object,here::here("data","Processed_Datasets","ledergorMyeloma","ledergorMyeloma_default","ledergorMyeloma_default.Rds"))
