slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
print(slurm_arrayid)

library(Seurat)
library(GloScope)
library(reticulate)
library(sceasy)
reticulate::use_condaenv("scvi-env", required=TRUE)

sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)

source("muscat_sim_serials.R")
#source("True_KLv2.R")
load("../../results/simulation/data/muscat_1sample_f5s.Rda")
x@rowRanges@elementMetadata$beta$cluster_id  =
  DataFrame(apply(x@rowRanges@elementMetadata$beta$cluster_id , 2, function(x) x/2))

set.seed(as.numeric(slurm_arrayid))
nk = 8
ns = 10
ng = 1e4
dd = TRUE
force = TRUE
rel_lfc = rep(1,nk)
paired = TRUE
lfc = 0.05
nc = 5000

probs = list(c(0.4, 0.05, 0.2, 0.05, 0.05, 0.05, 0.1, 0.1),
             c(0.15, 0.05, 0.1, 0.4, 0.05, 0.1, 0.05, 0.1),
             NULL,
             NULL)
p_dd = c(0.9, 0, 0.1, 0, 0, 0)



par_list = Step1(x, mysd = 0.2, nk,ns,ng, dd, force, rel_lfc, paired, lfc, probs, p_dd,  p_type = 0,
                 phylo_tree =NULL, offset = NULL, cats = cats_vec,  samp_off = NULL)
set.seed(as.numeric(slurm_arrayid))
rm(x)
sim = Step2(par_list = par_list, nc = nc, probs = probs, cats = cats_vec,p_ep = 0.5, p_dp = 0.3, p_dm = 0.5,
            off_var = NULL)


sim = as.Seurat(sim, data = NULL)
sim = RenameAssays(object = sim, originalexp = "RNA")

adata <- convertFormat(sim, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
scvi$data$setup_anndata(adata)
model_batch = scvi$model$SCVI(adata)

# train the model
model_batch$train()
latent_batch = model_batch$get_latent_representation()

# put it back in our original Seurat object
latent_batch <- as.matrix(latent_batch)
rownames(latent_batch) = colnames(sim)
colnames(latent_batch) <- paste0("scvi_", 1:ncol(latent_batch))


sim  = sim %>%
  NormalizeData(verbose = T, normalization.method = "LogNormalize")%>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)%>%
  ScaleData(verbose = T)%>%
  RunPCA(verbose = T)


plot_df = cbind(sim@meta.data, sim@reductions$pca@cell.embeddings, latent_batch)

set.seed(1)

dist_mat_GMM_pca = gloscope(x = plot_df, sample_id = "sample_id", dim_redu = "PC", ndim = 10, dens = "GMM",
                           BPPARAM = BiocParallel::MulticoreParam(2,RNGseed = 1), dist_mat = "KL", varapp = FALSE,
                           returndens = FALSE, epapp = FALSE)

set.seed(1)

dist_mat_KNN_pca = gloscope(x = plot_df, sample_id = "sample_id", dim_redu = "PC", ndim = 10, dens = "KNN",
                            BPPARAM = BiocParallel::MulticoreParam(2,RNGseed = 1), dist_mat = "KL", varapp = FALSE,
                            returndens = FALSE, epapp = FALSE)


set.seed(1)
dist_mat_GMM_scvi= gloscope(x = plot_df, sample_id = "sample_id", dim_redu = "scvi", ndim = 10, dens = "GMM",
                           BPPARAM = BiocParallel::MulticoreParam(2,RNGseed = 1), dist_mat = "KL", varapp = FALSE,
                           returndens = FALSE, epapp = FALSE)


set.seed(1)

dist_mat_KNN_scvi = gloscope(x = plot_df, sample_id = "sample_id", dim_redu = "scvi", ndim = 10, dens = "KNN",
                            BPPARAM = BiocParallel::MulticoreParam(2,RNGseed = 1), dist_mat = "KL", varapp = FALSE,
                            returndens = FALSE, epapp = FALSE)


save(dist_mat_KNN_pca,dist_mat_GMM_pca,
  dist_mat_KNN_scvi,dist_mat_GMM_scvi,
    file = paste0("../../results/simulation/distmat/baseline/sd020_lfc005/",slurm_arrayid,".Rda"))

