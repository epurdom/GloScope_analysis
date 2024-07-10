library(harmony)
library(SingleCellExperiment)
library(scuttle)
sce <- readRDS("../../data/Processed_Datasets/perezLupus/perezLupus_default/perezLupus_default.Rds")
library(scater)
myPCA = reducedDim(sce)
colnames(myPCA) = paste0("PC_",1:ncol(myPCA))
meta = as.data.frame(colData(sce))
meta$sample = as.character(meta$sample)

harmonized_pcs <- HarmonyMatrix(
  data_mat  = myPCA,
  meta_data = meta,
  vars_use  = "sample",
  do_pca    = FALSE
)

colnames(harmonized_pcs) = paste0("harmony_",1:ncol(harmonized_pcs))
save(harmonized_pcs, file = "../../results/Perez2022/harmonySam_pipeline.Rda")
