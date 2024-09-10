# GWAS meta-analysis of aPTT as a way to understand coagulation.
This repository contains supplementary material (scripts, figures, and tables) used during the realization of my Final Master Project, titled: GWAS Meta-analysis of aPTT as a Way to Understand Coagulation.

Workflow Diagram of the procedures conducted in this project: 

![GWAS Data Preprocessing (6)](https://github.com/user-attachments/assets/944da14b-6111-419e-8a9f-339f2594d491)

## Contents
- **Scripts**: Code used for data processing, statistical analysis, and visualization.
- **Figures**: Graphs and plots generated during the analysis.
- **Tables**: Data tables and summary statistics.

## Objective
The main goal of this project is to identify genetic variants associated with aPTT through a comprehensive meta-analysis of genome-wide association studies (GWAS). 

## Usage
The scripts and data provided here are intended for reproducibility and further analysis. Feel free to explore and adapt the material for related studies or projects.

# GWAS data preprocessing

The code used during the data harmonisation of the provided cohorts is included in the R folder with the name GWAS_data_processing.R

# Quality Control performed with EasyQC and Meta-analysis performed with METAL software

We first performed quality control using EasyQC package from R. Then, we performed the meta-analysis using METAL software.

Scripts used for this tools are located in the tools folder with the name EasyQC_script.ecf and MetalScript.txt

# Top tables

To process meta-analysis dataframes and extract the top SNPs, use the script located in the R folder: script_toptable_metaanalysis.R

# TWAS 

For TWAS analysis, follow these steps:

- Data Harmonization: Harmonize data from meta-analysis using the script located in the BASH folder: harmonization_script.sh
- Data Imputation: Impute data with the following script, also found in the BASH folder: imputation_script.sh
- Harmonized and Imputed Data Merging: Merge harmonized data with imputed data using the script: run_imputation.sh

Prediction and Meta-analysis (PrediXcan and MetaXcan softwares):
- prediXcan Software: Use the script spredixcan_script.R located in the R folder.
- MetaXcan Software: Use the script smultixcan_script.sh located in the BASH folder.

# Visualization
Visualization scripts are available in the plots folder:

PrediXcan z-score plot: zscoreplot_predixcan.R
Top SNPs Manhattan plot: topsnps_manhattan_plot.R


