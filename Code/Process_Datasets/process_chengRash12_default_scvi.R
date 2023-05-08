source("Code/Process_Datasets/run_scvi.R")

# run scvi for the default rash12 processing
# reduction keys must be set manually when updating the Seurat object
# see `run_scvi.R` for details about how to name reduction keys
seurat_object <- readRDS("data/Processed_Datasets/chengRash12/chengRash12_default/chengRash12_default.Rds") 

default_scvi <- run_scvi(seurat_object,batch_var=NULL)
seurat_object[["scvi"]] <- CreateDimReducObject(embeddings = default_scvi, key = "scvi_", assay = DefaultAssay(seurat_object))

saveRDS(seurat_object,"data/Processed_Datasets/chengRash12/chengRash12_default/chengRash12_default.Rds")
