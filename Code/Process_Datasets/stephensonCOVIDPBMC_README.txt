This folder contains the default Seurat processing of peripheral blood mononuclear cells (PBMCs) from COVID-infected patients and healthy controls presented in "Single-cell multi-omics analysis of the immune response in COVID-19" by Stephenson et al. (2021) [doi: 10.1038/s41591-021-01329-2]. 

The study contains 143 samples from 130 patients. There are a number of phenotypes including five infection severity levels and three types of control. Outside a handful of control units, the replicated samples are longitudinal at one of two additional stages of infection. The data was collected from three sample sites, which causes a marked batch effect.

An `.h5ad` file containing the processed counts in Seurat format was downloaded from the European Bioinformatics Institute with Accession E-MTAB-10026. The provided object includes PCA and Harmony embeddings we use as input to our distance calculations. We do not apply our usual quality control step of removing genes expressed in less than 10 cells. The object contains 647,366 cells and 24,929 genes.

Processed by: William Torous, Purdom Group
Date Downloaded: March 22nd, 2021
Date Processed: April 30th, 2022
