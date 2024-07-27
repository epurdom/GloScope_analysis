library(Seurat)
library(harmony)

seurat_object <- readRDS("../../data/Processed_Datasets/stephensonCOVIDPBMC/stephensonCOVIDPBMC_default/stephensonCOVIDPBMC_default.Rds")

V = seurat_object@reductions$pca@cell.embeddings

meta_data <- myseurat@meta.data
harmony_embeddings <- harmony::HarmonyMatrix(
  V, meta_data, 'Site', do_pca = FALSE, verbose=FALSE
)
plot_df = cbind(myseurat@meta.data,harmony_embeddings)
save(plot_df, file = "../../results/BatchStudy/rerun_harmonyallsample.Rda")