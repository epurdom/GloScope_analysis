This folder contains the default processing of Lung fibrosis from 6 studies. The associated paper is "Identification of a broadly fibrogenic macrophage subset induced by type 3 inflammation" by Fabre et al. (2023) [doi: 10.1126/sciimmunol.add8945].

The study design contains 144 samples from 6 different studies. For sick individuals, their phenotype is either chronic obstructive pulmonary disease, pulmonary fibrosis, or pulmonary systemic sclerosis.


UMI counts were downloaded from the SingleCell portal. The downloaded object has 714,923 cells expressing 19,982 genes. We note that this number of cells matches the post-quality-control total presented in the normalized data. Therefore, we do not apply further quality control beyond removing genes which are excluded from the normalized data; this leaves 19,680 genes. We apply our standard pre-processing steps for Seurat objects.

Processed by: Hao Wang, Purdom Group
Date Downloaded: April 11th, 2024
Date Processed: April 11th, 2024
