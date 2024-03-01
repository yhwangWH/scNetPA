Our code is divided into three parts:

1. normalizeAndGetIni.R is used for data normalization using Seurat, and the data generated is used as initial values for perturbation experiments.

2. SCENICForGRNInferrence.py is used for inferring Gene Regulatory Networks (GRN) using the SCENIC method. Some datasets required for preparation are included in the .py file itself, or you can refer to https://resources.aertslab.org/cistarget/ or the GitHub homepage at https://github.com/aertslab/pySCENIC.

2. scNetPA.R contains our scNetPA method, which includes network validation, network propagation, and perturbation analysis. The Protein-Protein Interaction (PPI) dataset used inside can be downloaded from the STRING dataset: https://figshare.com/articles/dataset/STRING_database_ver_10_5_archive_9606_data_sets/7857185.

The single-cell datasets for K562 and hepatocellular carcinoma used in the code can be downloaded from https://www.ncbi.nlm.nih.gov/geo/.
