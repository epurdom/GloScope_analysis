This folder contains the default Seurat processing of mouse brain cells presented in "A taxonomy of transcriptomic cell types across the isocortex and hippocampal formation" by Yao et al. (2021) [doi:10.1016/j.cell.2021.04.021]. A SingleCellExperiment with the quality-controlled UMI counts is downloaded using the R package `AllenInstituteBrainData`. The data set we consider is "Allen_Mouse_2020."

The study design "[profiles] ∼1.3 million cells covering the entire adult mouse isocortex and HPF and derived a transcriptomic cell-type taxonomy revealing a comprehensive repertoire of glutamatergic and GABAergic neuron types." The study design collected 59 samples from 15 brain regions across 54 mice. After filtering genes expressed in less than 10 cells, this yields 26,877 genes measured across 1,169,213 cells. Our standard pre-processing pipeline for SingleCellExperiments is used for processing, and this data set is too large to run scVI on.

Processed by: William Torous, Purdom Group                                   
Date Downloaded: June 14th, 2022
Date Processed: April 30th, 2023
