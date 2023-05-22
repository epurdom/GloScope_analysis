This folder contains the default processing of peripheral blood mononuclear cells (PBMCs) from 162 patients with systemic lupus erythematosus (SLE) and 99 healthy controls. The associated paper is "Single-cell RNA-seq reveals cell type–specific molecular and genetic associations to lupus" by Perez et al. (2022) [doi: 10.1126/science.abf1970].

The study design contains 274 samples from 261 patients. For sick individuals, their phenotype is either managed, flared, or post-flare treated. There are additionally 4 processing batches, and some samples were processed multiple times. We consider the cross of biological sample and processing cohort as units for sample-level analysis, and there are 336 of these.

UMI counts were downloaded from the CellxGene portal. The downloaded object has 1,263,676 cells expressing 30,933 genes. We note that this number of cells matches the post-quality-control total presented in Figure S1-C of the original paper. Therefore, we do not apply further quality control beyond removing genes which are expressed in less than 10 cells; this leaves 22,444 genes. We apply our standard pre-processing steps for SingleCellExperiment objects.

Processed by: William Torous and Hao Wang, Purdom Group
Date Downloaded: April 2nd, 2023
Date Processed: April 30th, 2023
