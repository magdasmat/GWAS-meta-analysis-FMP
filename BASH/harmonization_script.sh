#!/bin/bash

python path/to/gwas_parsing.py \
-gwas_file path/to/input_file.txt.gz \
-separator '\t' \
-output_column_map rsid variant_id \
-output_column_map chr chromosome \
-output_column_map pos position \
--chromosome_format \
--enforce_numeric_columns \
-output_column_map Allele2 non_effect_allele \
-output_column_map Allele1 effect_allele \
-output_column_map Effect effect_size \
-output_column_map Pvalue pvalue \
-output_column_map StdErr standard_error \
-output_column_map N sample_size \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size \
-snp_reference_metadata path/to/reference_metadata/gtex_v8_eur_filtered_maf0.01_monoallelic_variants_updated.txt.gz METADATA \
-output path/to/output_file.txt.gz
