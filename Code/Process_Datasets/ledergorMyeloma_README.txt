This folder contains the default Seurat processing of blood plasma cells from patients with multiple myeloma and healthy controls, presented in "Single cell dissection of plasma cell heterogeneity in symptomatic and asymptomatic myeloma" by Ledergor et al. (2018) [doi: 10.1038/s41591-018-0269-2]

The UMI counts were downloaded from the GEO using ascension number GSE117156. We remove two samples with less than 50 cells and one sample with extreme outlying distances after running our methodology. This leaves a design with 42 patients contributing 132 samples. There are individuals with 4 distinct types of multiple myeloma and healthy control blood plasma. Some of the sick patients provided samples before and after an oncology treatment. Two types of tumor samples are considered: circulating plasma cells (circPC) and bone marrow plasma cells (BMPC).

Following the original authors, we discard cells with less than 500 UMIs. We also remove genes expressed in fewer than 10 cells and run the default Seurat pipeline. This creates a final object consisting of 42,342 cells expressing 26,622 genes.

Processed by: William Torous, Purdom Group
Date Downloaded: November 10, 2022
Date Processed: May 2nd, 2023
