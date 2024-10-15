# DGRPool
Scripts used for generating figures of the DGRPool manuscript

# Citation
Please cite as: eLife2024;12:RP88981 DOI: https://doi.org/10.7554/eLife.88981.3

# Online Website
[https://dgrpool.epfl.ch/](https://dgrpool.epfl.ch/)

# Reproducibility
## Phenotyping data
In order to be fully reproducible, we downloaded the phenotypes on the website at a given timepoint. The script used, [1_download_phenotypes.R](1_download_phenotypes.R), access our API to download a "studies.json" file containing all metadata for each study. Then it uses the same API to download the phenotypes study by study, and format everything in a common format.
It then generates a RDS file with a given timestamp, which is the common file used by all other methods, so that all scripts are using the same data, collected at a given timestamp. We here provided the data used in the latest version of the manuscript in [data.all_pheno_15_07_24_filtered.rds](RDS/data.all_pheno_15_07_24_filtered.rds), but you can run again [1_download_phenotypes.R](1_download_phenotypes.R) to generate a new RDS with the latest up-to-date phenotyping data.

## Correlation
We also created a script, [2_correlation_phenotypes.R](2_correlation_phenotypes.R), to generate all correlations across phenotypes. This is generating correlation files in the [resources](resources) folder, which are then used by other scripts (such as Figure 2A, 2B, S5).

# Figures
All scripts needed to reproduce Figures of the paper are in their respective FigureXX folders. In general they are all using the timestamped RDS file that is present in the RDS folder (see Reproducibility section).

# GWAS
In the GWAS folder, we provide the two R scripts that are used to generate the GWAS results. Both for the pre-calculated GWAS ran for all phenotypes in the DGRPool database ([running_GWAS.R](GWAS/running_GWAS.R)), and for the user-submitted GWAS ([running_GWAS_user.R](GWAS/running_GWAS_user.R)). It requires [PLINK2](https://www.cog-genomics.org/plink/2.0/) to be installed on the machine to run correctly.
[Variant annotation](GWAS/dgrp.fb557.annot.txt.gz) and [covariate](GWAS/dgrp.cov.tsv) files are also in the GWAS folder.

# DGRP / Bloomington mapping
Since it can be hard to find DGRP/bloomington mapping files, especially for old DGRP lines, we provide a mapping file here: [resources/DGRP_Bloomington_Mapping.txt](resources/DGRP_Bloomington_Mapping.txt)

