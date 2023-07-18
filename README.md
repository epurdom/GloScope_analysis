# Replication Code for the GloScope Paper

This repository contains the code used to generate the figures, simulation experiments, and real-world data applications presented in the following paper:

[H. Wang, W. Torous, B. Gong, E. Purdom (2023).
Visualizing scRNA-Seq Data at Population Scale with GloScope. bioRxiv.](https://doi.org/10.1101/2023.05.29.542786)

The R code which generates the paper's figures can be found in `Code/Figure`. Likewise, the simulations can be rerun with the code in `Code/simulation`.
The folder `Code/Process_Datasets` contains the steps to transform publicly available scRNA-Seq UMI counts into a `Seurat` or `SingleCellExperiment` object with a standardized format.
Due to size constraints, the datasets used as input are not included in this repository. There are `.txt` files in that folder which link to the source of each dataset.  
With `Code/Analyze_Datasets`, the GloScope methodology is applied to these standardized objects, and the results are compared across datasets.

For further information about the R implementation of the GloScope methodology, visit [this GitHub repository]([url](https://github.com/epurdom/GloScope)https://github.com/epurdom/GloScope)https://github.com/epurdom/GloScope)https://github.com/epurdom/GloScope).
