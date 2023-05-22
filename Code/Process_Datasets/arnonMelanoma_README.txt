This folder contains the default Seurat processing of melanoma cells presented in "A Cancer Cell Program Promotes T Cell Exclusion and Resistance to Checkpoint Blockade" by Jerby-Arnon et al. (2018) [doi: 10.1016/j.cell.2018.09.006].

The dataset includes 7,186 malignant, stroma, and immune cells obtained from human melanoma tumors across two batches. Of interest are the characteristics of tumors unresponsive to immune checkpoint inhibitor (ICI) pharmacotherapy. The phenotypes are untreated or post-ICI (resistant); patient Mel04.3 is the only patient and sample who is post-ICI (responsive).

The already processed data was found via the GEO with accession number GSE115978. We note by inspection that all samples have reasonable ranges for cell counts, features per cell, and counts per cell. This reinforces that the data has indeed been quality controlled. So too does Table S2 in the supplementary materials which states there were originally 10,123 cells but 2,937 were low quality. Furthermore, no MT genes are present, suggesting they have been filtered already. The only additional QC step we take is to remove genes which appear in less than 10 cells. This leaves 21,493 features.

Supplementary Table S1 provides clarification that of the 33 samples, only 4 come from the 2 patients with replication (Mel129 and Mel75). The table also clarifies that the `orig.ident` covariate in the raw Seurat object is formatted with multiple duplicate patient encodings. We rectify this in our metadata processing. There are additional sample covariates in Table S2 we do not incorporate, such as location in the body and whether the tumor is primary or metastasized. The data set has two waves, an earlier pilot study (Tirosh et al., 2016) and the expansion presented in this paper; we consider this as a batch effect alongside patient source. 

The original authors uses Seurat to pre-process their data. Based on the `ImmRes_OE.R` file in the paper's GitHub repository, we infer that they apply the default pipeline parameters, which we do too. We also note that `ImmRes_Rfiles.zip` associated with this experiment in the Single Cell Portal contains a further QC data set with 6,173 cells and 10,483 genes. 

Processed by: William Torous, Purdom Group
Date Downloaded: November 20th, 2020
Date Processed: April 29th, 2023
