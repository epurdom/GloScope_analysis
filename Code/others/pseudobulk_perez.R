library(stringr)
library(Seurat)
library(SingleCellExperiment)
library(muscat)

sce = readRDS("../../data/Processed_Datasets/perezLupus/perezLupus_default/perezLupus_default.Rds")
sce$sample_id = paste0(sce$sample_uuid,"_", sce$Processing_Cohort)


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

plot_df <- cbind( seurat_sam@reductions$pca@cell.embeddings)


save(plot_df, file = "../../results/pseudobulk/pseudobulk_perez.Rda")

