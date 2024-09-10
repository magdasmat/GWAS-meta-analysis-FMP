#!/bin/bash

python  /path/to/gwas_summary_imputation_collect.py  \
-gwas_file /path/to/gwas_file_harmonized.txt.gz \
-folder /path/to/imputation_folder \
-pattern "chr*" \
-parsimony 7 \
-output /path/to/output/harmonized_imputed.txt.gz