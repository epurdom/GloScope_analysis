import scanpy as sc
import numpy as np
import scvi
import pandas as pd
adata = scvi.data.read_h5ad("../../data/COVID_143/covid_portal_210320_with_raw.h5ad")

scvi.data.setup_anndata(adata, layer="raw", batch_key="sample_id")
model = scvi.model.SCVI(adata)
model.train()
latent = model_scvi.get_latent_representation()

#np.savetxt('../../results/COVID_143/scvi_data.dat', latent)
pd.DataFrame(latent).to_csv("../../results/COVID143/scvi_sam.csv")

