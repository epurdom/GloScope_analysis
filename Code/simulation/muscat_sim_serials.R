library(muscat)
library(SingleCellExperiment)
library(varhandle)
library(dplyr)
library(purrr)
library(Seurat)
library(DirichletReg)
cats_vec = muscat:::cats # load DE category factors

#' Helper function to sample from a negative binomial distribution
#' 
#' @param cs value
#' @param d value
#' @param m value
#' @param lfc value
#' @param f value
#' @return 
#' @examples 
#' 
.nb <- function(cs, d, m, lfc = NULL, f = 1) {
  n_gs <- length(d)
  n_cs <- length(cs)
  if (is.null(lfc)) {
    lfc <- rep(0, n_gs)
  } else {
    lfc[lfc < 0] <- 0
  }
  fc <- f * (2 ^ lfc)
  fc_d <- fc
  fc <- rep(fc, each = n_cs)
  ds <- rep(1/d, each = n_cs)
  if(n_cs == 0){
    ms = NULL
    y = NULL
  }else{
#    ms <- c(t(matrix(replicate(n_cs,m),nrow = n_gs))) * fc 
    ms <- c(t(m)) * fc 
    y <- rnbinom(n_gs * n_cs, size = ds, mu = ms)
    y <- matrix(y, byrow = TRUE, 
                nrow = n_gs, ncol = n_cs, 
                dimnames = list(names(d), cs))
    ms <- split(ms, rep(seq_len(nrow(m)), each = n_cs))
    ms <- do.call("rbind", ms)
  }
  
  #  fc = split(fc, rep(seq_len(nrow(m)), each = n_cs))
  #  list(counts = y, means = ms, foldc = fc, origmean = m)
  list(counts = y, foldc = fc_d, means = ms)
}

#' 
#' 
#' @param cat value
#' @param "ep" value
#' @param "de" value
#' @param "dp" value
#' @param "dm" value
#' @param "db") value
#' @param cs_g1 value
#' @param cs_g2 value
#' @param m_g1 value
#' @param m_g2 value
#' @param d value
#' @param lfc value
#' @param ep value
#' @param dp value
#' @param dm value
#' @return 
#' @examples 
#' 
.sim <- function(
  cat = c("ee", "ep", "de", "dp", "dm", "db"),
  cs_g1, cs_g2, m_g1, m_g2, d, lfc, ep, dp, dm) {
  
  ng1 <- length(cs_g1)
  ng2 <- length(cs_g2)
  
  re <- switch(cat,
               "ee" = {
                 list(
                   .nb(cs_g1, d, m_g1),
                   .nb(cs_g2, d, m_g2))
               },
               "ep" = {
                 g1_hi <- sample(ng1, round(ng1 * ep))
                 g2_hi <- sample(ng2, round(ng2 * ep))
                 list(
                   .nb(cs_g1[-g1_hi], d, m_g1),
                   .nb(cs_g1[ g1_hi], d, m_g1, lfc), # 50% g1 hi
                   .nb(cs_g2[-g2_hi], d, m_g2),
                   .nb(cs_g2[ g2_hi], d, m_g2, lfc)) # 50% g2 hi
               },
               "de" = {
                 list(
                   .nb(cs_g1, d, m_g1, -lfc), # lfc < 0 => all g1 hi
                   .nb(cs_g2, d, m_g2,  lfc)) # lfc > 0 => all g2 hi
               },
               "dp" = {
                 props <- sample(c(dp, 1 - dp), 2)
                 g1_hi <- sample(ng1, round(ng1 * props[1]))
                 g2_hi <- sample(ng2, round(ng2 * props[2]))
                 list(                           
                   .nb(cs_g1[-g1_hi], d, m_g1), 
                   .nb(cs_g1[ g1_hi], d, m_g1,  lfc), # lfc > 0 => dp/(1-dp)% up
                   .nb(cs_g2[-g2_hi], d, m_g2), 
                   .nb(cs_g2[ g2_hi], d, m_g2, -lfc)) # lfc < 0 => (1-dp)/dp% up
               },
               "dm" = {
                 g1_hi <- sample(ng1, round(ng1 * dm))
                 g2_hi <- sample(ng2, round(ng2 * dm))
                 list(
                   .nb(cs_g1[-g1_hi], d, m_g1),
                   .nb(cs_g1[ g1_hi], d, m_g1, -lfc), # lfc < 0 => 50% g1 hi
                   .nb(cs_g2[-g2_hi], d, m_g2),
                   .nb(cs_g2[ g2_hi], d, m_g2,  lfc)) # lfc > 0 => 50% g2 hi
               }, 
               "db" = {
                 if (sample(c(TRUE, FALSE), 1)) {
                   # all g1 mi, 50% g2 hi
                   g2_hi <- sample(ng2, round(ng2 * 0.5))
                   list(
                     .nb(cs_g1, d, m_g1, abs(lfc), 0.5),
                     .nb(cs_g2[-g2_hi], d, m_g2, -lfc), 
                     .nb(cs_g2[ g2_hi], d, m_g2,  lfc)) 
                 } else {
                   # all g2 mi, 50% g1 hi
                   g1_hi <- sample(ng1, round(ng1 * 0.5))
                   list(
                     .nb(cs_g2, d, m_g2, abs(lfc), 0.5), 
                     .nb(cs_g1[-g1_hi], d, m_g1, -lfc),  
                     .nb(cs_g1[ g1_hi], d, m_g1,  lfc))  
                 }
               }
  )
  cs <- map(re, "counts")
  cs <- do.call("cbind", cs)
  ms <- map(re, "means") 
#      lapply(function(x) do.call(rbind, x))
  ms <- do.call("cbind", ms)
  colnames(ms) <- colnames(cs)
  rownames(ms) <- rownames(cs)
    rmv <- vapply(ms, is.null, logical(1))
  fc <- map(re, "foldc") 
  #  names(fc) = rownames(cs)
  #    lapply(function(x) do.call(rbind, x))
  #  fc = do.call("cbind", fc)
  #  colnames(fc) = colnames(cs)
  #  rownames(fc) = rownames(cs)
  
  #  o_m = map(re, "origmean")
  #  o_m = do.call("cbind", o_m)
  #  ms <- map(re, "means")
  #  names(ms) <- c("A", "B")[c(ng1, ng2) != 0]
  #  names(fc) <- c("A", "B")[c(ng1, ng2) != 0]
  fc <- do.call("cbind", fc)
  #  ms <- do.call("cbind", ms)
  list(cs = cs, ms = ms, fc = fc)
}


#' 
#' 
#' @param n value
#' @param ids value
#' @param probs value
#' @param pi_i value
#' @return 
#' @examples 
#' 
.sample_cell_md_new <- function(n, ids, probs = NULL, pi_i) {
  ns <- vapply(ids, length, numeric(1))
  if (is.null(probs)) 
    probs <- vector("list", 4)
  probs <- lapply(seq_along(probs), function(i) {
    if (!is.null(probs[[i]])) {
      return(probs[[i]])
    } else {
      rep(1 / ns[i], ns[i]*2)
    }
  })
  new_s <- paste0(rep(ids[[3]],2), rep(c(".A", ".B"), each = length(ids[[3]])))
  cd <- data.frame("cluster_id" = NA,
                   "sample_id" = sample(new_s, n, TRUE, probs[[3]]))
  cd$group_id <- substr(cd$sample_id, nchar(cd$sample_id), nchar(cd$sample_id)+1)
  #  probs_k <- rbind(rdirichlet(length(ids[[3]]), probs[[1]]*100),
  #                  rdirichlet(length(ids[[3]]), probs[[2]]*100))
  probs_k <- pi_i
#  rownames(probs_k) <- new_s
#  colnames(probs_k) <- ids[[1]]
  for(i in new_s){
    cd$cluster_id[cd$sample_id == i] <- sample(colnames(probs_k),
                                               sum(cd$sample_id == i),
                                               T,
                                               probs_k[i,])
  }
  
  cd$group_id <- factor(cd$group_id)
  cd$sample_id <- substr(cd$sample_id, 1,nchar(cd$sample_id)-2)
  return(cd)
}

#' Generate simulation parameters from fit SincleCellExperiment
#' 
#' This function uses the distributions fit to an input SincleCellExperiment
#' object to generate parameters for a specified simulation.
#' Only designs balanced on phenotype are supported.
#' 
#' @param x (SingleCellExperiment) Output of muscat::prepSim() with NB parameter estimates
#' @param mysd (numeric) Standard deviation of patient-level variability
#' @param nk (integer) The number of cell-type classes
#' @param ns (integer) The number of reference samples available
#' @param ng (integer) The number of genes
#' @param dd (logical) Whether to simulate differential expression
#' @param force (logical) Whether to force the specified # of samples
#' @param rel_lfc (numeric) A vector of scalings for the phenotypic fold-change by cluster. An atomic input is applied to all clusters
#' @param paired (logical) Whether or not to have a paired design
#' @param lfc (numeric) The lfc for DE genes in each cluster
#' @param probs (nested vector of numeric)
#' @param p_dd (vector of numeric) Unnormalized multinomial parameters for gene DE type
#' @param p_type (numeric) prob. of EE/EP gene being of class "type", c.f. "state" and equal genes in `muscat`
#' @param phylo_tree (logical) Whether to fit a hierarchical clustering during assignment of DE genes to clusters
#' @param offset (numeric or NULL) If numeric, the grand mean for offset in the NB model
#' @param cats (factor) The differential distribution category trypes
#' @param samp_off (numeric or NULL) If numeric, parameter for symmmetric uniform additive noise in offset
#' @param alpha (numeric) sample level cluster proportion variation
#' @return par_list (vector) Simulation parameters for Step2 of pipeline
#' @examples 
#' 
Step1 = function(x, mysd, nk, ns, ng, dd, force, rel_lfc,
	paired, lfc, probs, p_dd, p_type = 0,
	phylo_tree =NULL, offset = NULL,
	cats = cats_vec, samp_off = 0.5, alpha = 100){
	# Get unique clusters and sample IDs from SCE, as well as their count
	# kids0 and sids0 are named vectors with the value and name extracted
	nk0 <- length(kids0 <- rlang::set_names(levels(x$cluster_id)))
	ns0 <- length(sids0 <- rlang::set_names(levels(x$sample_id)))
	# a `muscat` fn. determines the number of simulated samples to match study design
	ns <- muscat:::.get_ns(ns0, ns, dd, paired, force)

	# Define standard alphanumeric IDs for sim (not SCE!) clusters, samples, and phenotypes
	nk <- length(kids <- rlang::set_names(paste0("cluster", seq_len(nk))))
	sids <- rlang::set_names(paste0("sample", seq_len(ns)))
	gids <- rlang::set_names(c("A", "B")) # two phenotypes

	# Map SCE group IDs to simulation slots defined immediately prior
	ref_kids <- rlang::set_names(sample(kids0, nk, nk > nk0), kids)
	if (!dd || paired) {
		# use same set of reference samples for both groups
		ref_sids <- sample(sids0, ns, ns > ns0)
		ref_sids <- cbind(ref_sids, ref_sids)
	} else {
		# draw reference samples at random for each group
		sidsA <- sample(sids0, ns, force && ns > ns0)
		if (force) {
			sidsB <- sample(sids0, ns, ns > ns0)
		} else {
			sidsB <- setdiff(sids0, sidsA)
			sidsB <- sample(sidsB, ns)
		}
		ref_sids <- cbind(sidsA, sidsB)
	}
	dimnames(ref_sids) <- list(sids, gids)

	# Generate list of relative logFCs for each cluster
	if (is.null(rel_lfc)) 
	rel_lfc <- rep(1, nk)
	if (is.null(names(rel_lfc))) {
		names(rel_lfc) <- kids
	} else {
		stopifnot(names(rel_lfc) %in% kids0)
	}

	# Draw samples' cluster frequencies using phenotype-specific Dirichlet parameters
	if(is.na(alpha)){
	  pi_i = rbind(matrix(probs[[1]], nrow = length(sids),ncol = nk, byrow = T),
	               matrix(probs[[2]], nrow = length(sids), ncol = nk,byrow = T))
	}else{
	  pi_i <- rbind(DirichletReg::rdirichlet(length(sids), probs[[1]]*alpha),
	                DirichletReg::rdirichlet(length(sids), probs[[2]]*alpha))
	}
	rownames(pi_i) <- paste0(rep(sids,2), rep(c(".A", ".B"), each = length(sids)))
	colnames(pi_i) <- kids

	# sample # of genes to simulate per category & gene indices
	n_dd <- table(sample(cats, ng, TRUE, p_dd))
	n_dd <- replicate(nk, n_dd)
	colnames(n_dd) <- kids
	gs <- paste0("gene", seq_len(ng))
	#While the number of EE and DE genes is the same between clusters, the exact genes vary
	gs_idx <- muscat:::.sample_gene_inds(gs, n_dd)

	# for each cluster, sample set of genes to simulate from
	gs_by_k <- rlang::set_names(sample(rownames(x), ng, ng > nrow(x)), gs)
	gs_by_k <- replicate(nk, gs_by_k)
	colnames(gs_by_k) <- kids

	# when 'phylo_tree' is specified, induce hierarchical cluster structure in `gs_by_k`
	if (!is.null(phylo_tree)) {
		res <- muscat:::.impute_shared_type_genes(x, gs_by_k, gs_idx, phylo_tree, phylo_pars)
		gs_by_k <- res$gs_by_k
		class  <- res$class
		specs <- res$specs
		# otherwise, simply impute type-genes w/o phylogeny
	} else if (p_type != 0) {
		res <- muscat:::.impute_type_genes(x, gs_by_k, gs_idx, p_type)
		stopifnot(!any(res$class == "shared"))
		gs_by_k <- res$gs_by_k
		class <- res$class
		specs <- res$specs
	} else {
		class <- rep("state", ng)
		specs <- rep(NA, ng)
		names(class) <- names(specs) <- gs
	}

	# Create a list of lists containing the name and number of genes in a cluster
	gs_by_k <- split(gs_by_k, col(gs_by_k))
	gs_by_k <- setNames(map(gs_by_k, set_names, gs), kids)
	gs_by_kc <- lapply(kids, function(k) 
		lapply(unfactor(cats), function(c) 
			gs_by_k[[k]][gs_idx[[c, k]]])) 

	# Create a nested matrix with the sampled log fold-change for each gene.
	# Fold-changes are sampled from Gamma with shape = 4 and rate = 4/lfc.
	# Signs are added with equal probability.
	lfc <- vapply(kids, function(k) 
		lapply(unfactor(cats), function(c) { 
			n <- n_dd[c, k]
			if (c == "ee") return(rep(NA, n))
			signs <- sample(c(-1, 1), n, TRUE)
			lfcs <- rgamma(n, 4, 4/lfc) * signs
			names(lfcs) <- gs_by_kc[[k]][[c]]
			lfcs * rel_lfc[k]
		}), vector("list", length(cats)))
	rownames(lfc) = unfactor(cats)

	# Create a  dataframe with dummy encoding of NB beta additive coefficients
	bs <- rowData(x)$beta
	if (!is.null(bs$sample_id)) {
		bs$sample_id <- cbind(0, bs$sample_id)
		names(bs$sample_id) <- sids0
	}
	if (!is.null(bs$cluster_id)) {
		bs$cluster_id <- cbind(0, bs$cluster_id)
		names(bs$cluster_id) <- kids0
	}
	sids_new = paste0(rep(sids, each = 2),c(".A",".B"))
	# add patient level variability through Gaussian noise in beta
	bs$sample_id <- DataFrame(matrix(rnorm(length(bs$beta0)*ns*2, sd = mysd),ncol = ns*2))
	colnames(bs$sample_id) = sids_new

	# Generate NB offsets for each sample
	os <- x$offset
	if(is.null(offset)){
		os_mean = mean(os)
	}else{
		os_mean = offset
	}
	# Add uniform noise for variability in offset
	if(is.null(samp_off)){
		os_mean <- rep(os_mean, ns*2)
	}else{
		os_mean <- os_mean+ runif(ns*2, -samp_off, samp_off)
	}
	names(os_mean) <- paste0(rep(sids,2), rep(c(".A", ".B"), each = length(sids)))

	b0 <- bs$beta0
	ds <- rowData(x)$disp
	rowData(x)$beta$sample_id = bs$sample_id

	# split input cells by cluster-sample
	cs_by_ks <- muscat:::.split_cells(x)

	# create a long data.frame for each gene,cluster cross and its NB parameters
	gi <- data.frame(
		gene = unlist(gs_idx),
		cluster_id = rep.int(rep(kids, each = length(cats)), c(n_dd)),
		category = rep.int(rep(cats, nk), c(n_dd)),
		logFC = unlist(lfc),
		sim_gene = unlist(gs_by_kc),
		sim_disp = ds[unlist(gs_by_kc)]) %>%
			mutate_at("gene", as.character)

	# sort genes alphabetically
	o <- order(as.numeric(gsub("[a-z]", "", gi$gene)))
	gi <- gi[o, ]; rownames(gi) <- NULL

	# return the parameter list
	par_list =  list(pi = pi_i,
		bs = bs,
		refk = ref_kids,
		refs = ref_sids,
		ds = ds,
		lfc = lfc,
		gs_by_kc = gs_by_kc,
		gs_idx = gs_idx,
		os_mean = os_mean,
		n_dd = n_dd,
		gi = gi,
		cs_by_ks = cs_by_ks)

	return(par_list)
}

#' Simulate data using Step1 parameter
#'
#' More details to come...
#' 
#' @param par_list (vector) Output of `Step1` with simulation parameters
#' @param nc (integer) Number of cells
#' @param probs value
#' @param cats value
#' @param p_ep value
#' @param p_dp value
#' @param p_dm value
#' @param off_var value
#' @return 
#' @examples 
#' 
Step2 = function(par_list, nc, probs, cats = cats_vec,
	p_ep = 0.5, p_dp = 0.3, p_dm = 0.5, off_var = 0.1){
	# Get number of samples to simulate
	ns <- dim(par_list$pi)[1]/2 # assumes a balanced design
	# Generate simulation IDs
	nk <- length(kids <- set_names(paste0("cluster", seq_len(dim(par_list$pi)[2]))))
	sids <- set_names(paste0("sample", seq_len(ns)))
	gids <- set_names(c("A", "B"))

	# sample reference clusters & samples
	ref_kids <- par_list$refk
	ref_sids <- par_list$refs

  ng = length(unique(par_list$gi$gene))
  nc = nc*ns*2
  # initialize count matrix
  gs <- paste0("gene", seq_len(ng))
  cs <- paste0("cell", seq_len(nc))
  y <- matrix(0, ng, nc, dimnames = list(gs, cs))

  
  # sample cell metadata
  cd <- .sample_cell_md_new(
    n = nc, probs = probs,
    ids = list(kids, kids,sids, gids),
    pi_i = par_list$pi)
  rownames(cd) <- cs
  cd$cell_info <- NULL
  
  cs_idx <- muscat:::.split_cells(cd, by = colnames(cd))
  n_cs <- modify_depth(cs_idx, -1, length)
  
  # split input cells by cluster-sample
  cs_by_ks <-  par_list$cs_by_ks
  
  # sample nb. of genes to simulate per category & gene indices
  gs_idx <- par_list$gs_idx
  
  # for ea. cluster, sample set of genes to simulate from
  gs_by_kc <- par_list$gs_by_kc


  # sample logFCs
  lfc <- par_list$lfc

    # compute NB parameters
  bs <- par_list$bs
  b0 <- bs$beta0
  sids_new = paste0(rep(sids, each = 2),c(".A",".B"))
  os_mean = par_list$os_mean
  ds <- par_list$ds
  n_dd <- par_list$n_dd
  for (k in kids) {
    for (s in sids) {
      # get output cell indices
      s_new = paste0(rep(s,each = 2), c(".A", ".B"))
      ci <- cs_idx[[k]][[s]]
      
      # get reference samples, clusters & cells
      s0 <- ref_sids[s, ]
      k0 <- ref_kids[k]
      cs0 <- cs_by_ks[[k0]][s0]
      
      
      # sample cells to simulate from
      if(is.null(n_cs[[k]][[s]])){
        n_cs[[k]][[s]][1:2] = 0
      }
      cs_g1 <- sample(cs0[[1]], n_cs[[k]][[s]][[1]], TRUE)
      cs_g2 <- sample(cs0[[2]], n_cs[[k]][[s]][[2]], TRUE)
      
      # get NB parameters
      bs_ks_A <- cbind(b0,
                       bs$cluster_id[[k0]],
                       bs$sample_id[[s_new[1]]])
      bs_ks_B <- cbind(b0,
                       bs$cluster_id[[k0]],
                       bs$sample_id[[s_new[2]]])
      if(is.null(off_var)){
        os_g1 <- os_mean[s_new[1]]  + rep(0, length(cs_g1))
        os_g2 <- os_mean[s_new[2]] + rep(0, length(cs_g2))
      }else{
        os_g1 <- os_mean[s_new[1]] + runif(length(cs_g1),-off_var, off_var)
        os_g2 <- os_mean[s_new[2]] + runif(length(cs_g2),-off_var, off_var)
      }
      cd[unlist(ci), "cell_info"] <- c(cs_g1, cs_g2)
      
      for (c in cats[n_dd[, k] != 0]) {
        
        # get reference genes & output gene indices
        gs0_sub <- gs_by_kc[[k]]
        names(gs0_sub) = unfactor(cats)
        gs0 <- unlist(gs0_sub[c])
        gi <- gs_idx[[c, k]]
        
        # get NB parameters
        ds_kc <- ds[gs0]
        lfc_kc <- lfc[[c, k]]
        bs_ksc_A <- exp(rowSums(bs_ks_A[gs0, , drop = FALSE]))
        bs_ksc_B <- exp(rowSums(bs_ks_B[gs0, , drop = FALSE]))

        ms_g1 <- outer(bs_ksc_A, exp(os_g1), "*")
        ms_g2 <- outer(bs_ksc_B, exp(os_g2), "*")
        
        re <- .sim(c, cs_g1, cs_g2, ms_g1, ms_g2, ds_kc, lfc_kc, p_ep, p_dp, p_dm)
        y[gi, unlist(ci)] <- re$cs

      }
    }
  }
  # gene info
  gi <- data.frame(
    gene = unlist(gs_idx),
    cluster_id = rep.int(rep(kids, each = length(cats)), c(n_dd)),
    category = rep.int(rep(cats, nk), c(n_dd)),
    logFC = unlist(lfc),
    sim_gene = unlist(gs_by_kc),
    sim_disp = ds[unlist(gs_by_kc)]) %>% 
    mutate_at("gene", as.character)
  
  
  # reorder
  o <- order(as.numeric(gsub("[a-z]", "", gi$gene)))
  gi <- gi[o, ]; rownames(gi) <- NULL
  
  # cell info
  cd$group_id <- droplevels(cd$group_id)
  cd$sample_id <- if (dd) {
    factor(paste(cd$sample_id, cd$group_id, sep = "."))
  } else factor(cd$sample_id)
  m <- match(levels(cd$sample_id), cd$sample_id)
  gids <- cd$group_id[m]

  o <- order(gids)
  sids <- levels(cd$sample_id)[o]
  cd <- cd %>% 
    mutate_at("sample_id", factor, levels = sids) %>% 
    mutate_at("cluster_id", factor, levels = kids)
  if (!dd) {
    cd$group_id <- NULL
    ref_sids <- ref_sids[, 1]
  }
  ei <- data.frame(sample_id = sids, group_id = gids[o])
  md <- list(
    experiment_info = ei,
    n_cells = table(cd$sample_id),
    gene_info = gi,
    ref_sids = ref_sids,
    ref_kids = ref_kids,
    ref_beta = par_list$s,
#    ref_off = par_list$offset,
    os_mean = os_mean,
    ref_prob = par_list$pi)
  #    fc_info = fc_matrix,
  #    cell_info = cell_info)
  SingleCellExperiment(
    assays = list(counts = as.matrix(y)),
    colData = cd, metadata = md)
}
