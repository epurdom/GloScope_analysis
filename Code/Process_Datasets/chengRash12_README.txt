This folder contains the default Seurat processing of skin cells presented in "Transcriptional Programming of Normal and Inflamed Human Epidermis at Single-Cell Resolution" by Cheng et al. (2018) [doi: 10.1016/j.celrep.2018.09.006].

The dataset includes 92,889 human epidermal cells from 9 normal and 3 inflamed skin samples. The original authors consider each of the 12 samples as a batch. Those from the trunk, scalp, or foreksin are healthy controls while those with psoriasis are inflamed.

The data was downloaded from the European Genome-Phenome Archive with associated Study ID EGAS00001002927. There are two datasets listed here: EGAD00001004367 with technology "Illumina NovaSeq 6000" and EGAD00010001620 with technology "single cell RNA-seq." Both datasets require permission from the UCSF Cheng/Cho Lab, which Hao received. Hao used the EGA Download Client V3 to download the data, which did not give an option to specify which EGAD to use. Therefore, we conclude there is only one copy of the dataset.

The number of cells in the raw reads matches the number reported after QC in the original paper's "Supplemental Information" section. Therefore, we do not apply any additional QC steps beyond removing genes expressed in less than 10 cells. This leaves 92,889 cells and 19,769 genes in the final object.
 
Processed by: William Torous, Purdom Group
Date Downloaded: July 10th, 2022
Date Processed: April 29th, 2023
