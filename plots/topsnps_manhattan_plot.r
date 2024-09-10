# top snps manhattan plot

library(data.table)
library(ggrepel)
library(ggplot2)
library(cowplot)
library(dplyr)
library(biomaRt)

# Load data
results <- fread("path/to/toptable.txt")

# Define known genes
known_genes <- data.frame(
  rsid = c("rs138260833", "rs698078", "rs2289252", "rs2066854", "rs1801020", 
           "rs9268951", "rs2726950", "rs505922", "rs927826", "rs7135039", 
           "rs11064497", "rs139097404"),
  gene_name = c("Z99572.1", "AC068631.1, KNG1", "F11, F11-AS1", "FGG", "F12, GRK6", 
                "HLA-DRB9", "SCARA5", "ABO", "PLXDC2", "VWF", "C1S", "CATSPER2")
)

# Rename columns
results <- results %>%
  rename(pvalue = Pvalue, chromosome = chr, start_location = pos)

# Create gene_name column
results$gene_name <- known_genes$gene_name[match(results$rsid, known_genes$rsid)]

# BioMart: Add gene information
addgene <- function(res) {
  ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  res$gene_id <- substr(res$gene_name, 1, 15)  # Assuming 'gene_name' is used; adjust if needed
  ref_gene <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name',
                                   'start_position', 'end_position'), mart = ensembl)
  ref <- ref_gene[ref_gene$ensembl_gene_id %in% res$gene_id, ]
  res <- res[res$gene_id %in% ref_gene$ensembl_gene_id, ]
  if (identical(res$gene_id, ref[na.omit(match(res$gene_id, ref$ensembl_gene_id)), ]$ensembl_gene_id)) {
    ref2 <- ref[na.omit(match(res$gene_id, ref$ensembl_gene_id)), ]
    res$chromosome <- ref2$chromosome_name
    res$start_location <- ref2$start_position
  }
  return(res)
}

results <- addgene(results)

# Data processing
results$pvalue <- as.numeric(results$pvalue)
results$log_pval <- -log10(results$pvalue)
results <- results %>% filter(!is.na(log_pval))

results$chromosome <- as.numeric(results$chromosome)
ylimit <- max(results$log_pval, na.rm = TRUE) + 1
no <- nrow(results)
Sig_Z_Thresh <- -log10(0.05 / no)

# Data for plot
d <- results[order(results$chromosome, results$start_location),]
d <- d[!is.na(d$log_pval),]
d$pos <- NA
ticks <- NULL
lastbase <- 0
numchroms <- length(unique(d$chromosome))
d$pos <- as.numeric(d$pos)
d$chromosome <- as.numeric(d$chromosome)
for (i in unique(d$chromosome)) {
  if (i == 1) {
    d[d$chromosome == i, ]$pos <- d[d$chromosome == i, ]$start_location
  } else {
    lastbase <- lastbase + tail(subset(d, chromosome == i-1)$start_location, 1)
    d[d$chromosome == i, ]$pos <- d[d$chromosome == i, ]$start_location + lastbase
  }
  ticks <- c(ticks, d[d$chromosome == i, ]$pos[floor(length(d[d$chromosome == i, ]$pos) / 2) + 1])
}
ticklim <- c(min(d$pos), max(d$pos))
mycols <- rep(c("gray35", "gray72"), 60)
d$Sig_Z_Thresh <- Sig_Z_Thresh
d_sig <- d[which(d$log_pval > d$Sig_Z_Thresh),]
d_sig <- d_sig[rev(order(d_sig$log_pval)),]
d_sig <- d_sig[!duplicated(d_sig$gene_name),]

if (sum(d_sig$log_pval > 0) > 0) {
  d_sig_pos <- d_sig[d_sig$log_pval > 0,]
}

chr_labs <- as.character(unique(d$chromosome))
chr_labs <- as.numeric(chr_labs)
chr_labs[chr_labs == 19 | chr_labs == 21] <- ' '

# Generate plot
if (dim(d_sig)[1] == 0) {
  p <- ggplot(d, aes(x = pos, y = log_pval, colour = factor(chromosome))) +
    geom_point(size = 0.5) +
    scale_x_continuous(name = "Chromosome", breaks = ticks, labels = chr_labs) +
    scale_y_continuous(name = '-log10(p-value)', limits = c(0, ylimit)) +
    scale_colour_manual(values = mycols, guide = FALSE) +
    geom_hline(yintercept = Sig_Z_Thresh, colour = "blue")
  
} else {
  p <- ggplot(d, aes(x = pos, y = log_pval, colour = factor(chromosome))) +
    geom_point(size = 0.5) +
    scale_x_continuous(name = "Chromosome", breaks = ticks, labels = chr_labs) +
    scale_y_continuous(name = '-log10(p-value)', limits = c(0, ylimit)) +
    scale_colour_manual(values = mycols, guide = FALSE) +
    geom_hline(yintercept = Sig_Z_Thresh, colour = "blue") +
    geom_point(data = d_sig, aes(x = pos, y = log_pval), colour = "red", fill = 'red', size = 1.5)
  
  if (sum(d_sig$log_pval > 0) > 0) {
    p <- p + geom_text_repel(data = d_sig_pos, aes(x = pos, y = log_pval, label = gene_name), colour = 'black', nudge_y = 1, size = 2.5, force = 5, segment.alpha = 0.25, ylim = c(Sig_Z_Thresh + 0.1, NA))
  }
}

p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               axis.line = element_line(colour = "black"),
               axis.text.x = element_text(angle = 45, size = 8, hjust = 1))

# Save plot as PNG
png(filename = "path/to/output/meta_ea_manhattan_pvalues.png", width = 1024, height = 750, res = 100)
print(p)
dev.off()
