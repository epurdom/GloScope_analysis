This folder contains the default Seurat processing of rash skin samples presented in "Classification of human chronic inflammatory skin disease based on single-cell immune profiling" by Liu et al. (2022)[DOI: 10.1126/sciimmunol.abl9165]. Data is downloaded from: https://zenodo.org/record/6471748#.ZFvfEHbMLEZ, Manuscript Object.rds.

The study contains immune cells from 31 rash samples, each from a distinct patient. The phenotypes breakdown as: 7 atopic dermatitis (AD), 8 psoriasis vulgaris (PV), 2 lichen planus (LP), 1 bullous pemphigoid (BP), 6 clinical/histopathologically indeterminate rashes, and 7 healthy controls. Combined they offer 145,810 cells and 22,253 genes. The UMI data is obtained as a Seurat object from the original authors. Since this has already been pre-processed, we only additionally remove genes expressed in less than 10 cells.

Processed by: William Torous, Purdom Group
Date Downloaded: July 8th, 2022
Date Processed: April 29th, 2023
