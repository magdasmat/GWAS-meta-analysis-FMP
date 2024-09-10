library(stringr)
library(data.table)
library(dplyr)
library(tidyr)

getwd()
setwd("path/to/working_directory")

# data <- fread("")
meta <- fread("path/to/meta_analysis.TBL")

head(meta)
# Separate columns and conserve original
data <- separate(data, col = cptid, into = c("chr", "pos", "rm", "rv"), sep = ":", remove = FALSE)
# remove "chr" 
data$chr <- gsub("chr","", data$chr)

# to remove chr "X" from autosomal GWAS files
# > any(meta$chr == "X")
# [1] TRUE
# meta <- subset(meta, chr != "X")

# arrange columns by chr and pos
meta <- arrange(meta, chr, as.numeric(pos))

class(meta$`P-value`) # if its as character i have to change pval column to numeric before filtering. 
# It's better to do it while pre-processing samples to be used for locuszoom and also to get the top table
meta$`P-value` <- as.numeric(meta$`P-value`)

nrow(meta[meta$`P-value`==0,]) #5, if there are pvalues=0 we can change them to the minimum value R accept which is 5e-324
meta[meta$`P-value`==0,]$`P-value`<- 5e-324



filtered_meta <- meta[meta$pvalue < 5e-9,]
nrow(filtered_meta)
table(filtered_meta$chr)
Pvalue
toptable <- c()
for (i in c(1:22)){
  RESULTSCHR <- filtered_meta[filtered_meta$chr == i, ]
  while(nrow(RESULTSCHR) > 0){
    minPv <- min(RESULTSCHR$pvalue)
    range <- c(RESULTSCHR[RESULTSCHR$pvalue == minPv , ]$pos - 1000000 ,RESULTSCHR[RESULTSCHR$pvalue == minPv , ]$pos + 1000000)
    if (length(range) > 2){ # In case the same p-value is shared between different variants.
      range <- range[c(1, length(range)/2+1)]
    } else {range <- range}
    toptable <- rbind(toptable, RESULTSCHR[RESULTSCHR$pvalue == minPv , ][1,])
    RESULTSCHR <- RESULTSCHR[RESULTSCHR$pos < range[1] | RESULTSCHR$pos > range[2],]
  }
}

write.table(toptable, "toptable_meta.txt", sep = "\t", row.names = FALSE, quote = FALSE)