slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
print(slurm_arrayid)

library(Seurat)
library(popPackage)
library(reticulate)
library(sceasy)
reticulate::use_condaenv("scvi-env", required=TRUE)

sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)

clus_KL <- function(prop1, prop2){
  KLdist <- 0
  for(i in 1:length(prop1)){
    KLdist <- KLdist + prop1[i]*(log(prop1[i]) - log(prop2[i]))
  }
  return(KLdist)
}


cluster_function = function(cluster_table){
  clusprop = matrix(cluster_table, ncol = ncol(cluster_table), dimnames = dimnames(cluster_table))
  clusprop[which(clusprop==0)] <- 0.5
  clusprop <- t(apply(clusprop, 1, function(x) x/sum(x)))
  
  sample_names <- rownames(clusprop)
  all_combn <- t(combn(sample_names, 2))
  dist_vec <- c()
  
  for (i in 1:nrow(all_combn)){
    s1 <- all_combn[i, 1]
    s2 <- all_combn[i, 2]
    KL <- clus_KL(prop1 = clusprop[s1,], prop2 = clusprop[s2,]) +
      clus_KL(prop1 = clusprop[s2,], prop2 = clusprop[s1,])
    dist_vec <- c(dist_vec, KL)
  }
  
  
  clusprop_diss <- matrix(0, ncol = length(sample_names), nrow = length(sample_names))
  
  rownames(clusprop_diss) <- sample_names
  colnames(clusprop_diss) <-  sample_names
  
  for (i in 1:nrow(all_combn)){
    clusprop_diss[all_combn[i, 1], all_combn[i, 2]] <- dist_vec[i]
    clusprop_diss[all_combn[i, 2], all_combn[i, 1]] <- dist_vec[i]
  }
  
  return(clusprop_diss)
  
}


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

A_clus = c(100,10,5,10,15,10,5,10)
B_clus = c(25,10,5,10,15,10,5,10)
A_clus = A_clus/sum(A_clus)
B_clus = B_clus/sum(B_clus)


probs = list(A_clus,
             B_clus,
             NULL,
             NULL)
p_dd = c(1, 0, 0, 0, 0, 0)



par_list = Step1(x, mysd = 0.13, nk,ns,ng, dd, force, rel_lfc, paired, lfc, probs, p_dd,  p_type = 0,
                 phylo_tree =NULL, offset = NULL, cats = cats_vec,  samp_off = NULL, alpha = 10)
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

mytable_true = table(sim$sample_id, sim$cluster_id)
mytable_est_pca = table(sim$sample_id, sim$seurat_clusters)

sim[["scvi"]] <- CreateDimReducObject(embeddings = latent_batch, key = "scvi_", assay = DefaultAssay(sim))

sim = sim %>%
  FindNeighbors(reduction = "scvi", verbose = T) %>%
  FindClusters(verbose = T) 
mytable_est_scvi = table(sim$sample_id, sim$seurat_clusters)


plot_df = cbind(sim@meta.data, sim@reductions$pca@cell.embeddings, latent_batch)

set.seed(1)

dist_mat_GMM_pca = distMat(x = plot_df, sample_id = "sample_id", dim_redu = "PC", ndim = 10, dens = "GMM",
                           BPPARAM = BiocParallel::MulticoreParam(2,RNGseed = 1), dist_mat = "KL", varapp = FALSE,
                           returndens = FALSE, epapp = FALSE)

set.seed(1)

dist_mat_KNN_pca = distMat(x = plot_df, sample_id = "sample_id", dim_redu = "PC", ndim = 10, dens = "KNN",
                            BPPARAM = BiocParallel::MulticoreParam(2,RNGseed = 1), dist_mat = "KL", varapp = FALSE,
                            returndens = FALSE, epapp = FALSE)


set.seed(1)
dist_mat_GMM_scvi= distMat(x = plot_df, sample_id = "sample_id", dim_redu = "scvi", ndim = 10, dens = "GMM",
                           BPPARAM = BiocParallel::MulticoreParam(2,RNGseed = 1), dist_mat = "KL", varapp = FALSE,
                           returndens = FALSE, epapp = FALSE)


set.seed(1)

dist_mat_KNN_scvi = distMat(x = plot_df, sample_id = "sample_id", dim_redu = "scvi", ndim = 10, dens = "KNN",
                            BPPARAM = BiocParallel::MulticoreParam(2,RNGseed = 1), dist_mat = "KL", varapp = FALSE,
                            returndens = FALSE, epapp = FALSE)

clus_true = cluster_function(mytable_true)
clus_seurat_pca = cluster_function(mytable_est_pca)
clus_seurat_scvi = cluster_function(mytable_est_scvi)

save(clus_true,clus_seurat_pca, clus_seurat_scvi, 
     file =  paste0("../../results/simulation/distmat/clusprop_V3/Setting1/25_100_alpha10/cluster/",slurm_arrayid ,".Rda" ))

save(dist_mat_KNN_pca,dist_mat_GMM_pca,
  dist_mat_KNN_scvi,dist_mat_GMM_scvi, 
    file = paste0("../../results/simulation/distmat/clusprop_V3/Setting1/25_100_alpha10/Glos/",slurm_arrayid,".Rda"))

