min_count = 1
min_cells = 10
min_genes = 100
min_size = 50
group_keep = NULL

myseurat = readRDS("../../data/Processed_Datasets/stephensonCOVIDPBMC/stephensonCOVIDPBMC_default/stephensonCOVIDPBMC_default.Rds")

sce <- as.SingleCellExperiment(myseurat, assay = "raw")
sce <- sce[,sce$group == "Covid"]
sce <- sce[,sce$batch == "Ncl"]



library(edgeR)
library(muscat)
library(SingleCellExperiment)


#pat = unique(sce$patient_id)
#sub_pat =  sample(pat,5)
#print(sub_pat)
# to ensure being able to replicate code, here are the 
# 5 random samples we picked:

sub_pat = c("MH8919327", "MH9143277", "MH9179824", "MH9143326", "MH9143271
")
sce = sce[, sce$patient_id %in% sub_pat]


colnames(colData(sce))[10] = "cluster_id"
#colnames(colData(sce))[25] = "sample_id"
colnames(colData(sce))[16] = "group_id"
sce$group_id = droplevels(sce$group_id)
sce$sample_id = droplevels(sce$sample_id)
sce$cluster_id = droplevels(sce$cluster_id)
x = sce
rm(sce)
colData(x) = colData(x)[, c("cluster_id", "sample_id", "group_id")]

############## below is following muscat's prepSim.R source code to 
############## preprocess and estimate parameters for simulation

vars <- c("sample_id", "cluster_id")
names(vars) <- vars <- intersect(vars, names(colData(x)))

# assure these are factors
for (v in vars) {
  # drop singular variables from model
  n <- length(unique(x[[v]]))
  if (n == 1) {
    rmv <- grep(v, vars)
    vars <- vars[-rmv]
  }
  if (!is.factor(x[[v]]))
    x[[v]] <- as.factor(x[[v]])
  x[[v]] <- droplevels(x[[v]])
}

n_cells0 <- ncol(x)
x <- muscat:::.update_sce(x)
if (is.null(group_keep)) {
  if ("group_id" %in% colnames(colData(x))) {
    group_keep <- levels(x$group_id)[1]
    cells_keep <- x$group_id == group_keep
  } else {
    cells_keep <- seq_len(ncol(x))
  }
} else {
  stopifnot(is.character(group_keep), 
            group_keep %in% levels(x$group_id))
  cells_keep <- x$group_id %in% group_keep
}
x <- x[, cells_keep]
x <- muscat:::.update_sce(x)

# keep genes w/ count > `min_count` in at least `min_cells`;
# keep cells w/ at least `min_genes` detected genes
genes_keep <- rowSums(counts(x) > min_count) >= min_cells
cells_keep <- colSums(counts(x) > 0) >= min_genes

x <- x[genes_keep, cells_keep, drop = FALSE]

# keep cluster-samples w/ at least 'min_size' cells
if (!is.null(min_size)) {
  n_cells <- table(x$cluster_id, x$sample_id)
  n_cells <-  muscat:::.filter_matrix(n_cells, n = min_size)
  if (ncol(n_cells) == 1)
    stop("Current 'min_size' retains only 1 sample,\nbut",
         " mean-dispersion estimation requires at least 2.")

  x <-  muscat:::.filter_sce(x, rownames(n_cells), colnames(n_cells))
}


f <- "~ 1"
for (v in vars)
  f <- paste(f, v, sep = "+")
cd <- as.data.frame(droplevels(colData(x)))
mm <- model.matrix(as.formula(f), data = cd)

# fit NB model

y <- DGEList(counts(x))
y <- calcNormFactors(y)
y <- estimateDisp(y, mm)
y <- glmFit(y, prior.count = 0)


cs <- y$coefficients
i <- !rowAnyNAs(cs)
cs <- cs[i, , drop = FALSE]
x <- x[i, , drop = FALSE]
ds <- y$dispersion[i]

# group betas by variable
bs <- DataFrame(
  beta0 = cs[, 1],
  row.names = rownames(x))
for (v in vars) {
  pat <- paste0("^", v)
  i <- grep(pat, colnames(cs))
  df <- DataFrame(cs[, i])
  nms <- colnames(cs)[i]
  names(df) <- gsub(pat, "", nms)
  bs[[v]] <- df
}
rowData(x)$beta <- bs

# store dispersions in row- & offsets in colData
names(ds) <- rownames(x)
rowData(x)$disp <- ds
os <- c(y$offset)
names(os) <- colnames(x)
x$offset <- os
save(x, file = "../../results/simulation/data/inter_muscat_5sample.Rda")

x@rowRanges@elementMetadata$beta$sample_id <- DataFrame(MH8919327 = 0)
x = x[,x$sample_id == "MH8919327"]

x$sample_id = droplevels(x$sample_id)
save(x, file = "../../results/simulation/data/muscat_1sample_f5s.Rda")

