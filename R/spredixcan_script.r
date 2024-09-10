library(data.table)
library(stringr)
library(magrittr)
library(glue)
library(plyr)

tissues <- c("Artery_Aorta", "Artery_Coronary", "Artery_Tibial",
             "Liver", "Whole_Blood")

launch_tissue <- function(tissue) {
    mashr_folder <- "/path/to/refpanels/TWAS_ref/ucsc/models/eqtl/mashr/"
    base_folder <- "/path/to/base_folder"
    
    orden <- glue("python /path/to/SPrediXcan.py ",
       "--gwas_file {base_folder}/DB/EA_harmonized_imputed.txt.gz ",
       "--model_db_path {mashr_folder}mashr_{tissue}.db ",
       "--covariance {mashr_folder}mashr_{tissue}.txt.gz ",
       "--snp_column panel_variant_id ",
       "--effect_allele_column effect_allele ",
       "--non_effect_allele_column non_effect_allele ",
       "--zscore_column zscore ",
       "--keep_non_rsid ",
       "--output_file {base_folder}/DB/spredixcan/EA_{tissue}.csv ",
       "--additional_output ",
       "--model_db_snp_key varID ",
       "--throw")
    
    # Execute order
    system(orden)
}

plyr::l_ply(tissues, launch_tissue)