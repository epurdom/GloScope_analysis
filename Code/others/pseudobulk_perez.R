library(stringr)
library(Seurat)
library(SingleCellExperiment)
library(muscat)

sce = readRDS("../../data/Processed_Datasets/perezLupus/perezLupus_default/initial_sce.Rds")
#sce$sample_id = paste0(sce$sample_uuid,"_", sce$Processing_Cohort)
sce@assays@data@listData$counts@seed@filepath = "../../data/Processed_Datasets/perezLupus/raw_files/scRNA_raw_cleaned.h5ad"
counts = counts(sce)
metadata = colData(sce)
sce = SingleCellExperiment(list(counts=counts), colData = metadata)

pb <- aggregateData(sce, assay = "counts", 
                    fun = "sum", 
                    by = c("sample"))
names(pb@assays) = "counts"

counts = pb@assays@data@listData$counts
seurat_sam <- CreateSeuratObject(counts)

seurat_sam  = seurat_sam %>%
  NormalizeData(verbose = T, normalization.method = "LogNormalize")%>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)%>%
  ScaleData(verbose = T)%>%
  RunPCA(npcs = min(c(ncol(seurat_sam)-1,50)),verbose = T)

plot_df <- cbind( seurat_sam@reductions$pca@cell.embeddings)


save(plot_df, file = "../../results/pseudobulk/pseudobulk_perez.Rda")

