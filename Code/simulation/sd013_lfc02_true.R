slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
print(slurm_arrayid)

load("../../results/COVID_143/muscat_1sample_f5s.Rda")
#library(popPackage)
source("muscat_sim_serials.R")
source("True_KL_serials_update.R")
x@rowRanges@elementMetadata$beta$cluster_id  =
  DataFrame(apply(x@rowRanges@elementMetadata$beta$cluster_id , 2, function(x) x/2))

set.seed(as.numeric(slurm_arrayid))
nk = 8
ns = 1
ng = 1e4
dd = TRUE
force = TRUE
rel_lfc = rep(1, nk)
paired = TRUE
lfc = 0.2
#nc = 5e4
#prob_s = runif(ns*2, 0.5, 1)
#prob_s = prob_s/sum(prob_s)

probs = list(c(0.4, 0.05, 0.2, 0.05, 0.05, 0.05, 0.1, 0.1),
             c(0.15, 0.05, 0.1, 0.4, 0.05, 0.1, 0.05, 0.1),
             NULL,
             NULL)
p_dd = c(0.9, 0, 0.1, 0, 0, 0)
mysd = 0.13

test = Step1(x, mysd = mysd,nk,ns,ng, dd, force, rel_lfc, paired, lfc, probs, p_dd,  p_type = 0, phylo_tree =NULL, offset = NULL, cats = cats, samp_off = NULL)
dist_mat = true_KL(par_list = test, nc = 10000, BPPARAM = BiocParallel::MulticoreParam(4))
save(dist_mat, file = paste0("../../results/simulation/vstrue/trueKL_Step1/sd013_lfc02/",slurm_arrayid ,".Rda"))

