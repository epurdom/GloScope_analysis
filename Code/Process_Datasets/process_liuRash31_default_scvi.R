source("Code/Process_Datasets/run_scvi.R")

# run scvi for the default rash31 processing
# reduction keys must be set manually when updating the Seurat object
# see `run_scvi.R` for details about how to name reduction keys
seurat_object <- readRDS("data/Processed_Datasets/liuRash31/liuRash31_default/liuRash31_default.Rds") 

default_scvi <- run_scvi(seurat_object,batch_var=NULL)
seurat_object[["scvi"]] <- CreateDimReducObject(embeddings = default_scvi, key = "scvi_", assay = DefaultAssay(seurat_object))

saveRDS(seurat_object,"data/Processed_Datasets/liuRash31/liuRash31_default/liuRash31_default.Rds")
