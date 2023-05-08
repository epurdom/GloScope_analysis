source("Code/Process_Datasets/run_scvi.R")

# run scvi for the default arnon processing
# reduction keys must be set manually when updating the Seurat object
# see `run_scvi.R` for details about how to name reduction keys
seurat_object <- readRDS("data/Processed_Datasets/arnonMelanoma/arnonMelanoma_default/arnonMelanoma_default.Rds") 

default_scvi <- run_scvi(seurat_object,batch_var=NULL)
seurat_object[["scvi"]] <- CreateDimReducObject(embeddings = default_scvi, key = "scvi_", assay = DefaultAssay(seurat_object))

batch_scvi <- run_scvi(seurat_object,batch_var="batch")
seurat_object[["scvi.batch"]] <- CreateDimReducObject(embeddings = batch_scvi, key = "scvi.batch_", assay = DefaultAssay(seurat_object))

patient_scvi <- run_scvi(seurat_object,batch_var="patient")
seurat_object[["scvi.pat"]] <- CreateDimReducObject(embeddings = patient_scvi, key = "scvi.pat_", assay = DefaultAssay(seurat_object))

saveRDS(seurat_object,"data/Processed_Datasets/arnonMelanoma/arnonMelanoma_default/arnonMelanoma_default.Rds")
