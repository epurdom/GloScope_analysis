This folder contains the default Seurat processing of lung cells presented in "A molecular single-cell lung atlas of lethal COVID-19" by Melms et al. (2021) [doi: 10.1038/s41586-021-03569-1].

The dataset contains "single-nucleus RNA sequencing of about 116,000 nuclei from the lungs of nineteen individuals who died with COVID-19 and underwent rapid autopsy and seven control individuals." One of the infected patients provided two samples, and all the other patients provided just one. The processed data was downloaded from the Broad Institute's Single Cell Portal. We confirm that the filtering steps described in the paper's Supplementary Materials match the provided data. We additionally remove genes expressed in less than 10 cells. The Supplementary Materials section also notes the default Seurat parameters were used before PCA, which we match. The final quality-controlled object contains 116,313 cells and 29,925 genes.

Processed by: William Torous, Purdom Group
Date Downloaded: July 8th, 2022
Date Processed: April 29th, 2023
