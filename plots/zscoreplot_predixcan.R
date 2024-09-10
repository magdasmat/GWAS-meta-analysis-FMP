library(data.table)
library(ggrepel)
library(ggplot2)
library(cowplot)
library(dplyr)
library(biomaRt)
options(ggrepel.max.overlaps=Inf)

# Create plotting
results1 <- fread("path/to/predixcan/output/file/Artery_Aorta.csv")
results2 <- fread("path/to/predixcan/output/file/Artery_Coronary.csv")
results3 <- fread("path/to/predixcan/output/file/Artery_Tibial.csv")
results4 <- fread("path/to/predixcan/output/file/Liver.csv")
results5 <- fread("path/to/predixcan/output/file/Whole_Blood.csv")

# BioMart:
#addgene <- function(res){
#ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
#res$gene_id=substr(res$gene,1,stop=15)
#ref_gene <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name',
#                               'start_position','end_position'), mart = ensembl)
#ref <- ref_gene[ref_gene$ensembl_gene_id %in% res$gene_id,]
#res <- res[res$gene_id %in% ref_gene$ensembl_gene_id,]
#if(identical(res$gene_id,ref[na.omit(match(res$gene_id,ref$ensembl_gene_id)),]$ensembl_gene_id)){
#  ref2 <- ref[na.omit(match(res$gene_id,ref$ensembl_gene_id)),];res$chromosome <- ref2$chromosome_name;res$start_location <- ref2$start_position}
#return(res)}
#results9 <- addgene(results9)
#fwrite(results9, "/home/gerard/CHARGE/TWAS/PSTALL/Input/mex-mashr-g1_Whole_Blood.csv", sep = ",")


#Merging MetaXcan results
merge=rbind(results1,results2,results3,results4,results5)
test=merge
test$var_g=NULL
test$largest_weight=NULL
library(dplyr)
test$z=abs(test$zscore)
test1=dplyr::group_by(test,gene_name) %>% summarise(max_zscore=max(z))
test2=dplyr::left_join(test1,test,by=c('gene_name'='gene_name','max_zscore'='z'))
test2=unique(test2)
test2$count=1
test3=dplyr::group_by(test2,gene_name,count) %>% summarise(count=sum(count))
test2$max_zscore=NULL

dataframe=test2
dataframe$chromosome <- as.numeric(dataframe$chromosome)
ylimit=max(abs(dataframe$zscore),na.rm=T)+1
no=nrow(test)
Sig_Z_Thresh=qnorm(1-(0.05/no)/2)
Sig_Z_Thresh=qnorm(1-(0.05/length(dataframe$zscore))/2)

d=dataframe[order(dataframe$chromosome, dataframe$start_location),]
d=d[!is.na(d$pvalue),]
d$pos=NA
ticks=NULL
lastbase=0
numchroms=length(unique(d$chromosome))
d$pos=as.numeric(d$pos)
d$chromosome <- as.numeric(d$chromosome)
for (i in unique(d$chromosome)) {
  if (i==1) {
    d[d$chromosome==i, ]$pos=d[d$chromosome==i, ]$start_location}	
  else {
    lastbase=lastbase+tail(subset(d,chromosome==i-1)$start_location, 1)
    d[d$chromosome==i, ]$pos=d[d$chromosome==i, ]$start_location+lastbase
  }
  ticks=c(ticks, d[d$chromosome==i, ]$pos[floor(length(d[d$chromosome==i, ]$pos)/2)+1])
}
ticklim=c(min(d$pos),max(d$pos))
mycols=rep(c("gray35","gray72"),60)
d$Sig_Z_Thresh<-Sig_Z_Thresh
d_sig<-d[which(abs(d$zscore) > d$Sig_Z_Thresh),]
d_sig<-d_sig[rev(order(abs(d_sig$zscore))),]
d_sig<-d_sig[!duplicated(d_sig$gene_name),]

if(sum(d_sig$zscore > 0) > 0){
  d_sig_pos<-d_sig[d_sig$zscore > 0,]
}
if(sum(d_sig$zscore < 0) > 0){
  d_sig_neg<-d_sig[d_sig$zscore < 0,]
}
chr_labs<-as.character(unique(d$chromosome))
chr_labs <- as.numeric(chr_labs)
chr_labs[chr_labs == '19'| chr_labs == '21']<-' '

if(dim(d_sig)[1] == 0){
  p<-ggplot(d,aes(x=pos,y=zscore,colour=factor(chromosome))) +
    geom_point(size=0.5) +
    scale_x_continuous(name="Chromosome", breaks=ticks, labels=chr_labs) +
    scale_y_continuous(name='Z score',limits=c(-ylimit,ylimit)) +
    scale_colour_manual(values=mycols, guide=FALSE) +
    geom_hline(yintercept=0,colour="black") +
    geom_hline(yintercept=Sig_Z_Thresh,colour="blue") +
    geom_hline(yintercept=-Sig_Z_Thresh,colour="blue")
     
} else {
  p<-ggplot(d,aes(x=pos,y=zscore,colour=factor(chromosome))) +
    geom_point(size=0.5) +
    scale_x_continuous(name="Chromosome", breaks=ticks, labels=chr_labs) +
    scale_y_continuous(name='Z score',limits=c(-ylimit,ylimit)) +
    scale_colour_manual(values=mycols, guide=FALSE) +
    geom_hline(yintercept=0,colour="black") +
    geom_hline(yintercept=Sig_Z_Thresh,colour="blue") +
    geom_hline(yintercept=-Sig_Z_Thresh,colour="blue") +
    geom_point(data=d_sig, aes(x=pos,y=zscore), colour="red", fill='red', size=1.5)
  
  if(sum(d_sig$zscore > 0) > 0){
    p<-p+geom_text_repel(data=d_sig_pos, aes(x=pos,y=zscore, label=gene_name), colour='black', nudge_y=1, size=2.5, force=5, segment.alpha=0.25, ylim=c(Sig_Z_Thresh+0.1,NA))
  }
          
  if(sum(d_sig$zscore < 0) > 0){
    p<-p+geom_text_repel(data=d_sig_neg, aes(x=pos,y=zscore, label=gene_name), colour='black', nudge_y=-1, size=2.5, force=5, segment.alpha=0.25, ylim=c(NA,-Sig_Z_Thresh-0.1))
  }
}
   
p<-p+theme(	panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(angle=45, size=8, hjust=1))

# Save as PNG
png(filename="path/to/output/file/predixcan_manhattanplot_pvalues.png",width=1024,height=750,res=100)
plot(p)
dev.off()



