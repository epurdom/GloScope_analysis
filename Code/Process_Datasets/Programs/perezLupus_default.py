import scanpy as sc
import pandas as pd
import anndata as ad
import os
import numpy as np
import anndata as ad

adata = sc.read_h5ad('../../data/Processed_Datasets/perezLupus/raw_files/local.h5ad')
#adata.layers["counts"] = adata_all.raw.X.copy()
#adata.layers["data"] = adata_all.X.copy()
adata.layers["counts"] = adata.raw.X.copy()
adata.layers["data"] = adata.X.copy()
adata.write('../../data/Processed_Datasets/perezLupus/raw_files/scRNA_raw_cleaned.h5ad')
