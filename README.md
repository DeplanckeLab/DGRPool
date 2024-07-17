# DGRPool
Scripts used for generating figures of the DGRPool manuscript

# Online Website
[https://dgrpool.epfl.ch/](https://dgrpool.epfl.ch/)

# Reproducibility
In order to be fully reproducible, we downloaded the phenotypes on the website at a given timepoint. The script used, [download_phenotypes.R](download_phenotypes.R), access our API to download a "studies.json" file containing all metadata for each study. Then it uses the same API to download the phenotypes study by study, and format everything in a common format.
It then generates a RDS file with a given timestamp, which is the common file used by all other methods, so that all scripts are using the same data, collected at a given timestamp. We here provided the data used in the latest version of the manuscript in [data.all_pheno_15_07_24_filtered.rds](RDS/data.all_pheno_15_07_24_filtered.rds), but you can run again [download_phenotypes.R](download_phenotypes.R) to generate a new RDS with the latest up-to-date phenotyping data.

# Figures
All scripts needed to reproduce Figures of the paper are in their respective FigureXX folders. In general they are all using the timestamped RDS file that is present in the RDS folder (see Reproducibility section).

# GWAS
In the GWAS folder, we provide the two R scripts that are used to generate the GWAS results. Both for the pre-calculated GWAS ran for all phenotypes in the DGRPool database ([running_GWAS.R](GWAS/running_GWAS.R)), and for the user-submitted GWAS ([running_GWAS_user.R](GWAS/running_GWAS_user.R)). It requires [PLINK2](https://www.cog-genomics.org/plink/2.0/) to be installed on the machine to run correctly.
[Variant annotation](GWAS/dgrp.fb557.annot.txt.gz) and [covariate](GWAS/dgrp.cov.tsv) files are also in the GWAS folder.

# DGRP / Bloomington mapping
Since it can be hard to find DGRP/bloomington mapping files, especially for old DGRP lines, we provide a mapping file here: [DGRP_Bloomington_Mapping.txt](DGRP_Bloomington_Mapping.txt)

