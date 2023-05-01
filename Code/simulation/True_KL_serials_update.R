library(muscat)
library(purrr)
library(varhandle)
library(SingleCellExperiment)
library(dplyr)
library(BiocParallel)
library(stringr)
#library(Seurat)
#library(DirichletReg)

#' 
#' 
#' @param nc value
#' @param kids value
#' @param gid value
#' @param sid value
#' @param pid value
#' @param bs value
#' @param lfc value
#' @param os_mean value
#' @param ds value
#' @return 
#' @examples 
#' 
sample_fun = function(nc, kids, gid, sid, pid, bs, lfc, os_mean, ds){
  cluster_id <- sample(kids, nc, TRUE, pid)

  y <- matrix(0, nrow(lfc), nc)
  bs_sample <- bs$sample_id[[sid]]
  for(i in kids){
    ci <- which(cluster_id == i)
    if(length(ci)>0){
      lfc_k <- lfc[[i]]
      lfc_k[is.na(lfc_k)] = 0
      if(gid == "A"){
        lfc_k = -lfc_k
      }
      lfc_k[lfc_k<0] = 0
      fc <- 1*(2 ^ lfc_k)
      bs_k <-  cbind(bs$beta0,
                     bs_sample,
                     bs$cluster_id[[i]])
      bs_k <- exp(rowSums(bs_k))
      ms <- outer(bs_k, exp(os_mean), "*")

      fc <- rep(fc, each =length(ci))
      d <- rep(1/ds, each = length(ci))
      m <- c(t(matrix(replicate(length(ci),ms),nrow = length(ms)))) * fc
      y[,ci] <- matrix(rnbinom(ng * length(ci), size = d, mu = m),byrow = TRUE,
                       nrow = ng, ncol = length(ci))
    }
  }
  return(y)
}

#' 
#' 
#' @param y value
#' @param kids value
#' @param bs value
#' @param sid value
#' @param gid value
#' @param lfc value
#' @param ds value
#' @param pid value
#' @param os_mean value
#' @return 
#' @examples 
#' 
int_dens <- function(y, kids, bs, sid, gid,lfc,ds, pid,os_mean){
  mydens = rep(NA, ncol(y))
  for(j in 1:ncol(y)){
    dens <- c()
    for(i in kids){
      bs_k <- cbind(bs$beta0,
                    bs$sample_id[[sid]],
                    bs$cluster_id[[i]])
      lfc_k <- lfc[[i]]
      lfc_k[is.na(lfc_k)] = 0
      if(gid == "A"){
        lfc_k <- -lfc_k
      }
      lfc_k[lfc_k<0] = 0
      fc <- 1*2^lfc_k
      bs_k <- exp(rowSums(bs_k))
      ms <- outer(bs_k, exp(os_mean), "*")
      m <- ms*fc
      dens <- c(dens, log(pid[i]) + sum(dnbinom(y[,j],size = 1/ds, mu = m, log = T)))
      }
    mydens[j] = max(dens) + log(sum(exp(dens - max(dens))))

    }

  return(mydens)
}

#' 
#' 
#' @param x value
#' @param nc value
#' @param sids value
#' @param kids value
#' @param gids value
#' @param pi_k value
#' @param bs value
#' @param lfc value
#' @param os_mean value
#' @param ds value
#' @return 
#' @examples 
#' 
calc_truedist = function(x, nc, sids, kids, gids, pi_k, bs, lfc, os_mean, ds){
  g_y <- gids[x]
  pi_y <- pi_k[x,]
  y <- sample_fun(nc, kids = kids, gid = g_y, sid = x, pid = pi_y, bs, lfc, os_mean, ds)
  dens_y <- int_dens(y, kids, bs, sid = x, gid = g_y, lfc, ds, pid = pi_y, os_mean)

  dist_y = c()
  for(j in sids){
    gid <- gids[j]
    pid <- pi_k[j,]
    dens_j <- int_dens(y, kids, bs, sid = j, gid = gid, lfc, ds, pid = pid, os_mean)
    dist_y = c(dist_y, mean(dens_y - dens_j))
  }
  return(dist_y)
}

#KL_fun <- function(nc, ep, kids, g1,g2, s1,s2, pi1,pi2, bs, lfc, os_mean, f){
#  y <- sample_fun(nc, ep, kids = kids, gid = g1, sid = s1, pi = pi1, bs, lfc, os_mean, f)
#  dens1 <- int_dens(y, kids, bs, sid = s1, gid = g1, lfc, ds, pi = pi1, os_mean, f)
#  dens2 <- int_dens(y, kids, bs, sid = s2, gid = g2, lfc, ds, pi = pi2, os_mean, f)

#  KL <- sum(log(dens1/(dens2+ep)))/length(dens1)
#  return(KL)
#}

#sKL_fun <- function(nc, ep, kids, g1,g2, s1,s2, pi1,pi2, bs, lfc, os_mean, f){
#  sKL <- KL_fun(nc, ep, kids, g1 = g1, g2 = g2, s1 = s1, s2 = s2, pi1 = pi1, pi2 = pi2, bs, lfc, os_mean, f) +
#    KL_fun(nc, ep, kids, g1 = g2, g2 = g1, s1 = s2, s2 = s1, pi1 = pi2, pi2 = pi1, bs, lfc, os_mean, f)
#  return(sKL)
#}



#' 
#' 
#' @param par_list value
#' @param nc value
#' @param BPPARAM value
#' @return 
#' @examples 
#' 
true_KL <- function(par_list, nc = 10000, BPPARAM = BiocParallel::MulticoreParam(4)){
  #get basic info
  ns <- dim(par_list$pi)[1]/2
  
  # simulation IDs
  nk <- length(kids <- set_names(paste0("cluster", seq_len(dim(par_list$pi)[2]))))
  sids <- colnames(par_list$bs$sample_id)
  gids <- str_sub(sids, nchar(sids), nchar(sids)+1)
  
#  sids <- sim@metadata$experiment_info$sample_id
#  gids <- sim@metadata$experiment_info$group_id
  names(gids) <- sids
  all_comb <- t(combn(sids, 2))
  ref_gene <- unique(par_list$gi[,c("gene", "sim_gene")])
  ng = length(unique(par_list$gi$gene))
  bs <- par_list$bs
  bs <- bs[ref_gene$sim_gene,]
  os_mean <- par_list$os_mean[1]
  b0 <- bs$beta0
  ds <- unique(par_list$gi[,c("gene", "sim_disp")])$sim_disp
  pi_k <- par_list$pi
  ref_kids <- par_list$refk
  colnames(bs$cluster_id) <- names(ref_kids)[match(colnames(bs$cluster_id),ref_kids)]

  KL_dist <- matrix(0, nrow = length(sids), ncol = length(sids))
  colnames(KL_dist) <- rownames(KL_dist) <- sids
#  g_info <- par_list$gi$logFC
  lfc <- reshape(par_list$gi[, c("gene", "cluster_id", "logFC")],
                 idvar = "gene", timevar = "cluster_id", direction = "wide")
#  lfc[is.na(lfc)] = 0
  lfc <- data.frame(lfc, row.names = 1)
  colnames(lfc) <- kids

  sids_list <- as.list(sids)
  dist = unlist(bplapply(sids_list, calc_truedist, nc=nc, sids=sids, kids=kids, gids=gids,
                pi_k=pi_k, bs=bs, lfc=lfc, os_mean=os_mean, ds=ds,
                BPPARAM = BPPARAM))


  dist_mat <- matrix(dist, nrow = length(sids), byrow = T)
  rownames(dist_mat) = colnames(dist_mat) = sids
  new_mat = matrix(NA, nrow = length(sids), ncol = length(sids))
  rownames(new_mat) = colnames(new_mat) = sids

  for(i in rownames(dist_mat)){
    for(j in colnames(dist_mat)){
      new_mat[i,j] = dist_mat[i,j] + dist_mat[j,i]
    }
  }
  diag(new_mat) = 0

  return(new_mat)
}
