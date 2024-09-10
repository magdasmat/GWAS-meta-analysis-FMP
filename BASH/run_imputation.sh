#!/bin/bash

# Define the list of chromosomes to process
chromosomes=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)

# Loop through each chromosome and execute the command
for chromosome in "${chromosomes[@]}"; do
  echo "Processing chromosome ${chromosome}..."
  python /path/to/gwas_summary_imputation.py \
    -by_region_file /path/to/by_region_file/TWAS_ref/eur_ld.bed.gz \
    -gwas_file /path/to/gwas_file_harmonized.txt.gz \
    -parquet_genotype_metadata /path/to/genotype_metadata/TWAS_ref/reference_panel_1000G/variant_metadata.parquet \
    -parquet_genotype /path/to/genotype/TWAS_ref/reference_panel_1000G/chr${chromosome}.variants.parquet \
    -window 100000 -parsimony 7 -chromosome ${chromosome} -regularization 0.1 -frequency_filter 0.01 \
    --use_palindromic_snps --standardise_dosages \
    -output /path/to/output_dir/chr${chromosome}.txt.gz
  echo "Chromosome ${chromosome} processed."
done
