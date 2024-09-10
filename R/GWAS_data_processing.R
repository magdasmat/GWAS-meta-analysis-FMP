# GWAS summary statistics data harmozation and inverse-normal transformation
# As each GWAs summary statistics has diferent units (some sec some rate). We will normalize them by dividing the betas by SD. 
# Here are three examples of the procedure followed to this aim. 
# load required packages

library(data.table)
library(dplyr)
library(tidyr)

# --------------------------------------------
# Japan GWAS Data Processing
# --------------------------------------------

# Load GWAS summary statistics for Japan
japan <- fread("path_to_GWAS_summary_statistics_file")  # GWAS data for aPTT from Japan Biobank
head(japan)
dim(japan)  # 5960951 rows
options(scipen = 999)  # Avoid scientific notation for better readability

# Create a new copy of the dataset
newjapan <- japan

# Add 'ID' column using the format "chr<CHR>"
newjapan$ID <- paste0("chr", japan$CHR)
head(newjapan)

# Add 'pos2' column (a copy of 'POS') and rearrange the columns
newjapan$pos2 <- japan$POS
japan2 <- newjapan %>% relocate(ID, .before = SNP)
japan2 <- japan2 %>% relocate(pos2, .before = SNP)
japan2 <- japan2 %>% relocate(POS, .before = pos2)
head(japan2)

# Save the processed data to a .bed file for LiftOver
write.table(japan2, "path_to_output_bed_file", append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Command to perform LiftOver to convert from hg19 to hg38
# Using the chain file: hg19 to hg38 conversion reference panel
# /path_to_liftOver -bedPlus=3 path_to_input_bed_file path_to_chain_file path_to_output_bed_file path_to_unlifted_file

# Load the LiftOver output (converted to hg38)
output <- fread("path_to_liftOver_output_file")
head(output)

# Rename columns to match the GWAS dataset
colnames(output) <- c("CHR", "pos", "post2", "rs", "chr", "ref", "alt", "EAF", "IMPUTATION", "BETA", "SE", "PVAL", "log10P", "N")
head(output)
dim(output)  # 5960951 rows

# Create 'cptid' column by combining 'CHR' and 'pos'
output$ID <- paste0(output$CHR, ":", output$pos)
names(output)[names(output) == 'ID'] <- 'cptid'
head(output)

# Rearrange columns to move 'cptid' before 'rs'
output <- output %>% relocate(cptid, .before = rs)

# Remove unnecessary columns ('CHR', 'pos', 'post2', 'chr')
japan2 <- output[, c(-1, -2, -3, -5)]
dim(japan2)  # 5960951 rows
head(japan2)
tail(japan2)

# Check for duplicated 'cptid'
sum(duplicated(japan2$cptid))  # 2 duplicated IDs

# Replace P-values of 0 with a small value to avoid errors in downstream analysis
japan3 <- japan2
japan2$PVAL <- ifelse(japan2$PVAL <= 0, 5e-323, japan2$PVAL)
table(identical(japan2$PVAL, japan3$PVAL))  # Check if replacement was successful
head(japan2)

# Save the processed data to a file
write.table(japan2, "path_to_final_output_file", append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

# Load the final dataset for further analysis
japan <- fread("path_to_final_output_file")
head(japan)

# Search for specific variants based on pattern
pattern <- "chr11:12541586"  # Example pattern to search for
matching_rows <- japan[grep(pattern, japan$cptid), ]

# Combine 'cptid' with effect and other allele information
japan$cptid_combined <- paste0(japan$cptid, ":", japan$EFFECT_ALLELE, ":", japan$OTHER_ALLELE)

# Rename columns to reflect that Japan uses the alternative allele as the effect allele
colnames(japan) <- c("cptid", "chr", "EFFECT_ALLELE", "OTHER_ALLELE", "EAF", "IMPUTATION", "BETA", "SE", "PVAL", "log10P", "N", "cptid_combined")
head(japan)

# Check for missing values
table(is.na(japan))

# Load reference allele frequency data (1000 Genomes ASN panel)
# Using the reference panel: 1000 Genomes ASN (Asian population)
a <- fread("path_to_reference_allele_frequency_file")
head(a)

# Sort both datasets by 'cptid'
ordered_data <- a[order(a$cptid), ]
ordered_japan <- japan[order(japan$cptid_combined), ]
head(ordered_japan)
head(ordered_data)

# Find common IDs between the datasets
common_ids <- intersect(ordered_japan$cptid_combined, ordered_data$cptid)
print(common_ids)
length(common_ids)  # 4129 common IDs

# Check dimensions of the datasets
dim(japan)
dim(a)

# Save the processed Japan dataset to a CSV file
write.csv(japan, "path_to_japan_csv_file", row.names = FALSE, sep = ",", quote = FALSE)

# Load necessary libraries for data processing
library(data.table)
library(dplyr)

# --------------------------------------------
# Middle Eastern GWAS Data Processing
# --------------------------------------------

# Load Middle Eastern GWAS summary statistics
middle <- fread("path_to_middle_eastern_GWAS_summary_statistics")
head(middle)
dim(middle)  # 8155707 rows

# Create a new dataset and adjust columns
newmiddle <- middle
newmiddle$ID <- paste0("chr", middle$chromosome)
newmiddle$pos2 <- middle$base_pair_location
newmiddle <- newmiddle %>% relocate(ID, .before = chromosome)
newmiddle <- newmiddle %>% relocate(base_pair_location, .after = ID)
newmiddle <- newmiddle %>% relocate(pos2, .after = base_pair_location)
head(newmiddle)

# Save the dataset in .bed format for LiftOver
write.table(newmiddle, "path_to_output_bed_file_for_middle", append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Run LiftOver to convert from hg19 to hg38 (requires .bed format)
# Reference panel: hg19 to hg38 chain file
# Command to run LiftOver:
# /path_to_liftOver -bedPlus=3 path_to_input_bed_for_middle path_to_chain_file path_to_output_bed path_to_unlifted_file

# Load the LiftOver output
output <- fread("path_to_liftOver_output_for_middle")
head(output)

# Rename columns to match the dataset structure
colnames(output) <- c("CHR", "pos", "post2", "chr", "p_value", "EFFECT_ALLELE", "OTHER_ALLELE", "EAF", "BETA")
head(output)
dim(output)  # 7879692 rows

# Prepare the final dataset
output$ID <- paste0(output$CHR, ":", output$pos)
names(output)[names(output) == 'ID'] <- 'cptid'
output <- output %>% relocate(cptid, .before = p_value)
output <- output[, c(-1, -2, -3, -4)]  # Remove unnecessary columns
output$N <- "5988"  # Sample size
output$SE <- sqrt((output$p_value * (1 - output$p_value)) / as.numeric(output$N))  # Calculate SE
head(output)

# Check for duplicated 'cptid'
sum(duplicated(output$cptid))  # 0 duplicates

# Save the processed data to CSV
write.csv(output, "path_to_final_output_file_for_middle", row.names = FALSE, sep = ",", quote = FALSE)

# Combine 'cptid' with alleles
output$cptid_combined <- paste0(output$cptid, ":", output$OTHER_ALLELE, ":", output$EFFECT_ALLELE)

# Example search pattern for specific variants
pattern <- "chr11:12541586"
matching_rows <- output[grep(pattern, output$cptid), ]
matching_rows

# ------------------------------------------------------------------------
# ARIC Cohort - European Ancestry (CHARGE cohorts, follow same procedure)
# ------------------------------------------------------------------------

# Load ARIC GWAS summary statistics
aric <- fread("path_to_ARIC_GWAS_summary_statistics")
head(aric)

# Adjust BETA and SE (rescale by lab reference factor 2.9)
aric$BETA <- aric$BETA / 2.9
aric$SE <- aric$SE / 2.9

# Create MAF column (minor allele frequency, optional)
aric$MAF <- ifelse(aric$AF_coded <= 0.5, aric$AF_coded, 1 - aric$AF_coded)
summary(aric$MAF)

# Filter data for MAF > 0.01
newaric <- aric[aric$MAF > 0.01, ]
dim(newaric)  # Filtered dataset, although this is performed by the quality control, so this is optional

# Save the processed ARIC dataset to CSV
write.csv(aric, "path_to_final_output_file_for_ARIC", row.names = FALSE, sep = ",", quote = FALSE)


