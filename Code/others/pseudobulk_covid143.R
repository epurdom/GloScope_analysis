library(stringr)
library(muscat)
library(Seurat)

myseurat = readRDS("../../data/Processed_Datasets/stephensonCOVIDPBMC/stephensonCOVIDPBMC_default/stephensonCOVIDPBMC_default.Rds")

sce <- as.SingleCellExperiment(myseurat, assay = "raw")
rm(myseurat)
sce <- prepSCE(sce, 
                   kid = "initial_clustering", # subpopulation assignments
                   gid = "Status",  # group IDs (ctrl/stim)
                   sid = "patient_id",   # sample IDs (ctrl/stim.1234)
                   drop = TRUE)


pb <- aggregateData(sce, assay = "counts", 
                           fun = "sum", 
                           by = c("sample_id"))
names(pb@assays) = "counts"

counts = pb@assays@data@listData$counts
seurat_sam <- CreateSeuratObject(counts)


seurat_sam  = seurat_sam %>%
  NormalizeData(verbose = T, normalization.method = "LogNormalize")%>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)%>%
  ScaleData(verbose = T)%>%
  RunPCA(npcs = min(c(ncol(seurat_sam)-1,50)),verbose = T)

plot_df <- cbind(seurat_sam@meta.data, seurat_sam@reductions$pca@cell.embeddings)


save(plot_df, file = "../../results/pseudobulk/pseudobulk_covid143.Rda")



