##### SF4 make PCOA and effect size plot for MLN RNA-seq samples


#### 7.16

library(reshape)
library(ggplot2)
library(matrixStats)

setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/data/MLN_RNA_Seq")

######
source("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/DIST_EffectSize/var_plotter.R")

rna_df <- read.table("normalizedCounts.txt", T)
meta <- read.table("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/data/metadata/mice_metadata.11.19_mouse_id.txt", 
                   header=T, sep='\t')
rownames(meta) <- as.character(meta$mouse_id)
colnames(rna_df) <- gsub("X", "", colnames(rna_df))
colnames(rna_df) <- gsub("_", "-", colnames(rna_df))
cyt_keep <- intersect(colnames(rna_df), meta$mouse_id)
meta <- meta[cyt_keep,]

input=log2(t(rna_df[,-1])+1)
colnames(input) <- rna_df[,1]
rownames(input) <- colnames(rna_df)[-1]

rowZeros <- apply(input, 1, function(x) {sum(x == 0)})
colZeros <- apply(input, 2, function(x) {sum(x == 0)})

#input = input[which(rowZeros<0.5*ncol(input)),which(colZeros<0.5*nrow(input)) ]
input = input[,which(colZeros<0.5*nrow(input)) ]

RV <- data.frame(vars=colVars(input), names=colnames(input))
top_var <- as.character(RV[order(RV$vars, decreasing=T),]$names[1:200])
input <- input[,top_var]

meta_keep<-meta[rownames(input),]

effectors=meta_keep[,c(2,3,4,6,8)]
col_var=meta_keep$Environment
shape_var=meta_keep$Genotype

gg_mln_rna <- var_plotter(input, effectors, col_var, "Environment",
                          shape_var, "Genotype", "MLN RNA-Seq")

pdf("mln_rnaseq_pcoa_ES.pdf", height = 5, width = 15)
grid.arrange(gg_mln_rna)
dev.off()


