#!/bin/bash

python /path/to/SMulTiXcan.py \
--models_folder /path/to/mashr_folder \
--models_name_pattern "mashr_(.*).db" \
--snp_covariance /path/to/ucsc/models/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz \
--metaxcan_folder /path/to/spredixcan \
--metaxcan_filter "EA_(.*).csv" \
--metaxcan_file_name_parse_pattern "EA_(.*).csv" \
--gwas_file /path/to/harmonized_imputed.txt.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--keep_non_rsid \
--model_db_snp_key varID \
--cutoff_condition_number 30 \
--verbosity 7 \
--throw \
--output /path/to/output/EA_smultixcan.txt
