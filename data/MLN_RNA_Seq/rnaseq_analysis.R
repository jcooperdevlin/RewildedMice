### RNASeq Analysis MLN
### 3.4

library(ggplot2)
library(DESeq2)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(scales)
library(ggrepel)
library(EnhancedVolcano)

#setwd("/Volumes/lokep01lab/lokep01labspace/Cooper_RNASeq/3.4")
setwd("/Volumes/research/lokep01lab/lokep01labspace/Rewilding_Data/int/data/MLN_RNA_Seq/")

meta <- read.table("metadata_mln.csv", T, ",")
meta$id <- gsub("_lab.*", "",meta$Sample.name)
meta$id <- gsub("_A.*", "",meta$id)
meta$id <- gsub("_N.*", "",meta$id)
meta$id <- gsub("_B.*", "",meta$id)
meta$id <- gsub("_", "-",meta$id)

orderer <-as.character(meta$id)

extra_meta <- read.table("mice_metadata.11.19_mouse_id.txt", T, '\t')
meta <- merge(meta, extra_meta, by.x = 'id', by.y='mouse_id', all.x=T)

write.table(meta, "metadata_full_RNASeq.txt", sep='\t', quote=F, row.names=F)

counts <- read.table("MLN_gene_exp.csv", T, ",")
count_fix <- gsub("_lab.*", "", colnames(counts))
count_fix <- gsub("_A.*", "", count_fix)
count_fix <- gsub("_N.*", "", count_fix)
count_fix <- gsub("_B.*", "", count_fix)
count_fix <- gsub("_", "-", count_fix)
count_fix <- gsub("X", "", count_fix)
colnames(counts) <- c("Gene", count_fix[-1])
rownames(counts) <- counts[,1]
counts <- counts[,-1]
counts <- counts[,as.character(meta$id)]
counts <- counts[1:(nrow(counts)-5),]

cds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ 1)
cds <- cds[ which( rowSums( counts( cds ) ) > 5 ), ]
dim( cds )

dds <- DESeq( cds )
norm_counts <- counts(dds, normalized = T)
norm_counts <- norm_counts[,orderer] ## for GEO
colnames(norm_counts) <- gsub("-", "_", colnames(norm_counts))
norm_counts <- data.frame(GeneSymbol=rownames(norm_counts), norm_counts)
colnames(norm_counts)[-1] <- gsub("X", "", colnames(norm_counts)[-1])


## organize for GEO
write.table(norm_counts, "normalizedCounts.txt", sep='\t', row.names=F, quote=F)

## all
norm_int <- norm_counts[order(rowVars(norm_counts),decreasing=T),][1:2000,]
pc <- prcomp(t(norm_int), scale.=T)
huh <- data.frame(pc$rotation)
pc_plot <- data.frame(meta, PC1 = pc$x[,1], PC2 = pc$x[,2])
dim(norm_int)

ggplot(pc_plot, aes(PC1,PC2, color = as.factor(Wedge_cage))) +
  geom_point(size=3)
pdf("Gender_pca.pdf", height = 4, width = 6)
ggplot(pc_plot, aes(PC1,PC2, color = Gender:Environment)) +
  geom_point(size=3)+ ggtitle(paste0(nrow(norm_int), " Variable Genes"))
dev.off()

g1=ggplot(pc_plot, aes(PC1,PC2, color = Environment, shape = Genotype)) +
  geom_point(size=3) + ggtitle(paste0(nrow(norm_int), " Variable Genes")) +
  scale_color_manual(values = c("red3", "mediumorchid3"))+theme_bw()
g2=g1+facet_wrap(~Genotype)
pdf("variableGenes_pca.pdf", height = 10, width = 8)
grid.arrange(g1, g1 + facet_wrap(~ Genotype, ncol = 2), ncol = 1)
dev.off()

data <- data.frame(obsnames=row.names(PC$x), PC$x[,1:2])
datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation[,1:2])
datapc$var1 <- rescale(datapc$PC1, c(min(data$PC1),max(data$PC1)))
datapc$var2 <- rescale(datapc$PC2, c(min(data$PC2),max(data$PC2)))

datapc$mult <- abs(datapc$PC1*datapc$PC2)
datapc <- datapc[order(datapc$mult, decreasing = T),]
datapc2 = datapc[1:25,]

g_seg1=ggplot(pc_plot, aes(PC1,PC2, color = Environment, shape=Genotype)) +
  geom_point(size=3) + ggtitle(paste0(nrow(norm_int), " Variable Genes")) +
  scale_color_manual(values = c("red3", "mediumorchid3")) + coord_equal()
g_seg=ggplot(pc_plot, aes(PC1,PC2)) +
  geom_point(size=3, color=NA) + ggtitle("Loadings for top 25 genes") +
  scale_color_manual(values = c("red3", "mediumorchid3")) + coord_equal() + 
  geom_segment(data=datapc2, aes(x=0, y=0, xend=var1, yend=var2), 
               arrow=arrow(length=unit(0.2,"cm")), alpha=0.5) +
  geom_text_repel(data=datapc2, aes(x=var1, y=var2, label=varnames), size = 3)
g_seg

pdf("variable_genes_PCA_biplot.pdf", height = 8, width = 10)
grid.arrange(g_seg1+theme(legend.position = 'none'), g_seg, nrow = 1)
dev.off()

## lab
lab_keep <- subset(meta, Environment == "lab")$id
lab_meta <- subset(meta, Environment == "lab")
lab_norm_counts <- norm_counts[,lab_keep]
norm_int <- lab_norm_counts[order(rowVars(lab_norm_counts),decreasing=T),][1:2000,]
pc <- prcomp(t(norm_int), scale.=T)
huh <- data.frame(pc$rotation)
pc_plot <- data.frame(lab_meta, PC1 = pc$x[,1], PC2 = pc$x[,2])
dim(norm_int)


g1=ggplot(pc_plot, aes(PC1,PC2, color = Environment, shape = Genotype)) +
  geom_point(size=3) + ggtitle(paste0(nrow(norm_int), " Variable Genes")) +
  scale_color_manual(values = c("red3", "mediumorchid3"))
g2=g1+facet_wrap(~Genotype)
pdf("Lab_variableGenes_PCA.pdf", height = 5, width = 10)
grid.arrange(g1+theme(legend.position = 'none'),g2,widths=c(1,1.3),ncol=2)
dev.off()

## wild
wild_keep <- subset(meta, Environment == "wild")$id
wild_meta <- subset(meta, Environment == "wild")
wild_norm_counts <- norm_counts[,wild_keep]
norm_int <- wild_norm_counts[order(rowVars(wild_norm_counts),decreasing=T),][1:2000,]
pc <- prcomp(t(norm_int), scale.=T)
huh <- data.frame(pc$rotation)
pc_plot <- data.frame(wild_meta, PC1 = pc$x[,1], PC2 = pc$x[,2])
dim(norm_int)


g1=ggplot(pc_plot, aes(PC1,PC2, color = Environment, shape = Genotype)) +
  geom_point(size=3) + ggtitle(paste0(nrow(norm_int), " Variable Genes")) +
  scale_color_manual(values = c("mediumorchid3"))
g2=g1+facet_wrap(~Genotype)
pdf("Wild_variableGenes_PCA.pdf", height = 5, width = 10)
grid.arrange(g1+theme(legend.position = 'none'),g2,widths=c(1,1.3),ncol=2)
dev.off()

#

#
#
#
#### DE Analysis
##
cds.subset <- DESeqDataSet( cds, design = ~ Environment )
dds.subset <- DESeq( cds.subset )
res <- data.frame(results( dds.subset, contrast = c("Environment", "wild", "lab")))
res$Gene <- rownames(res)
write.table(res, "Lab_v_wild_all.txt", quote = F, row.names=F, sep='\t')
res_filt <- subset(res, padj < 0.05)
res_filt <- subset(res, log2FoldChange < -1.5 | log2FoldChange > 1.5)

gene2keep <- as.character(res_filt$Gene)

#
#
#
### pause do DE by genotype

cds.subset <- DESeqDataSet( cds, design = ~ Genotype )
dds.subset <- DESeq( cds.subset )
contrasts_l <- list(
  c("Genotype", "AtgE", "AtgH"),
  c("Genotype", "AtgE", "B6"),
  c("Genotype", "AtgE", "NOD2"),
  c("Genotype", "AtgH", "B6"),
  c("Genotype", "AtgH", "NOD2"),
  c("Genotype", "B6", "NOD2")
)
res_all <- data.frame(res[1,])
res_all$comp <- NA
for(i in 1:6){
res <- data.frame(results( dds.subset, contrast = contrasts_l[[i]]))
res$Gene <- rownames(res)
res$comp <- paste0(contrasts_l[[i]][2], "_vs_", contrasts_l[[i]][3])
res_all <- rbind(res_all, res)
}
res_all <- res_all[-1,]
write.table(res_all, "Genotype_vs_all.txt", quote = F, row.names=F, sep='\t')

### for volcano colors

keyvals <- rep('grey50', nrow(res))
# set the base name/label as 'Mid'
names(keyvals) <- rep('NS', nrow(res))
# modify keyvals for transcripts with fold change > 2.5
keyvals[which(res$log2FoldChange > 1.2)] <- 'red3'
names(keyvals)[which(res$log2FoldChange > 1.2)] <- 'Wild'
# modify keyvals for transcripts with fold change < -2.5
keyvals[which(res$log2FoldChange < -1.2)] <- 'mediumorchid3'
names(keyvals)[which(res$log2FoldChange < -1.2)] <- 'NYU'
unique(names(keyvals))
##volcano
pdf("DE_genes_volcano.pdf", height = 7, width = 9)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                colAlpha =1,
                selectLab = rownames(res)[which(names(keyvals) %in% c('NYU', 'Wild'))],
                pCutoff = 0.05,
                FCcutoff = 1.2,
                transcriptLabSize = 4.0,
                colCustom = keyvals,
                xlim = c(-2, 6))
dev.off()

##


norm_counts <- counts(dds.subset, normalized = T)
int_genes <- unique(rownames(res_filt))
e_heat <- norm_counts[int_genes,]

va <- columnAnnotation(df = data.frame(Environment=meta$Environment,
                                       Genotype=meta$Genotype),
                    col = list(Environment=c(lab="red3",wild="mediumorchid3"),
                               Genotype=c(AtgE="green4", AtgH="dodgerblue2", B6="black", NOD2="orangered2")))

namer <- paste0("heats/DiffExp_Lab_Wild_", genotypes[i], ".pdf")
pdf("heats/DiffExp_Lab_Wild.pdf", height = 7, width = 5)
e_scale <- t(scale(t(e_heat)))
Heatmap(e_scale, top_annotation = va,
        heatmap_legend_param = list(title = "Scaled\nExpression"),
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        cluster_rows = T, cluster_columns = T,
        column_title = paste0(nrow(e_scale), " DiffExp Genes"),
        column_title_gp = gpar(fontsize = 15),
        show_column_names = F, show_row_names = T,
        row_names_gp = gpar(fontsize = 5))
dev.off()


res_geno <- data.frame(baseMean=NA,	log2FoldChange = NA,
                      lfcSE=NA,	stat=NA,	pvalue=NA,	padj=NA, Genes = NA, Genotype=NA)
genotypes <- as.character(unique(meta$Genotype))
for (i in 1:length(genotypes)){
  subber <- subset(meta, Genotype == genotypes[i])$id
  cds.subset <- cds[ , which( colnames( cds ) %in% subber ) ]
  meta_sub <- subset(meta, Genotype == genotypes[i])
  
  cds.subset <- DESeqDataSet( cds.subset, design = ~ Environment )
  dds.subset <- DESeq( cds.subset )
  res <- data.frame(results( dds.subset, contrast = c("Environment", "wild", "lab")))
  res$Genes <- rownames(res)
  res$Genotype <- rep(genotypes[i], nrow(res))
  res_geno <- rbind(res_geno, res)
  res_filt <- subset(res, padj < 0.05)
  res_filt <- subset(res, log2FoldChange < -2 | log2FoldChange > 2)
  write.table(res_filt, paste0("heats/Lab_v_Wild_", genotypes[i], "_DE.txt"), row.names=F, sep='\t', quote=F)
  
  norm_counts <- counts(dds.subset, normalized = T)
  int_genes <- unique(rownames(res_filt))
  e_heat <- norm_counts[int_genes,]
  
  va <- columnAnnotation(df = data.frame(Environment=meta_sub$Environment),
                         col = list(Environment=c(lab="red3",wild="mediumorchid3")))
  
  namer <- paste0("heats/DiffExp_Lab_Wild_", genotypes[i], ".pdf")
  pdf(namer, height = 7, width = 5)
  e_scale <- t(scale(t(e_heat)))
  pp=Heatmap(e_scale, top_annotation = va,
          heatmap_legend_param = list(title = "Scaled\nExpression"),
          col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
          cluster_rows = T, cluster_columns = T,
          column_title = paste0(nrow(e_scale), " DiffExp Genes"),
          column_title_gp = gpar(fontsize = 15),
          show_column_names = F, show_row_names = T,
          row_names_gp = gpar(fontsize = 5))
  print(pp)
  dev.off()
  
}
res_geno <- res_geno[-1,]
write.table(res_geno, "Lab_v_Wild_genotypes.txt", quote=F, row.names=F, sep='\t')





##
#
#
###### any genes overlapping USELESS*

# des <- list.files("heats", pattern = '.txt', full.names = T)
# 
# res_df <- data.frame(baseMean=NA,	log2FoldChange = NA,
#                        lfcSE=NA,	stat=NA,	pvalue=NA,	padj=NA, Genes = NA, Genotype=NA)
# for (i in 1:length(des)){
#   curr <- read.table(des[i], T, "\t")
#   res_df <- rbind(res_df, curr)
# }
# res_df <- res_df[-1,]
# 
# gene_tab <- data.frame(table(res_df$Genes))
# 
# res_merge <- merge(res_df, gene_tab, by.x = "Genes", by.y = "Var1", all.x = T)
# 
# res_int <- subset(res_merge, Freq > 1)
# good_genes <- unique(res_int$Genes)
# 
# 
# norm_counts <- counts(dds.subset, normalized = T)
# e_heat <- norm_counts[good_genes,]
# 
# va <- columnAnnotation(df = data.frame(Environment=meta$Environment,
#                                        Genotype=meta$Genotype),
#                        col = list(Environment=c(lab="red3",wild="mediumorchid3"),
#                                   Genotype=c(AtgE="green4", AtgH="dodgerblue2", B6="black", NOD2="orangered2")))
# 
# 
# pdf("heats/DiffExp_Lab_Wild_comboGenotypes.pdf", height = 7, width = 5)
# e_scale <- t(scale(t(e_heat)))
# Heatmap(e_scale, top_annotation = va,
#         heatmap_legend_param = list(title = "Scaled\nExpression"),
#         col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
#         cluster_rows = T, cluster_columns = T,
#         column_title = paste0(nrow(e_scale), " DiffExp Genes"),
#         column_title_gp = gpar(fontsize = 15),
#         show_column_names = F, show_row_names = T,
#         row_names_gp = gpar(fontsize = 5))
# dev.off()

PC=pc
data <- data.frame(obsnames=row.names(PC$x), PC$x[,1:2])
datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation[,1:2])
datapc$var1 <- rescale(datapc$PC1, c(min(data$PC1),max(data$PC1)))
datapc$var2 <- rescale(datapc$PC2, c(min(data$PC2),max(data$PC2)))

datapc$mult <- abs(datapc$PC1*datapc$PC2)
datapc <- datapc[order(datapc$mult, decreasing = T),]
datapc2 = datapc[datapc$varnames %in% gene2keep,]

g_seg1=ggplot(pc_plot, aes(PC1,PC2, color = Environment, shape=Genotype)) +
  geom_point(size=3) + ggtitle(paste0(nrow(norm_int), " Variable Genes")) +
  scale_color_manual(values = c("red3", "mediumorchid3")) + coord_equal()
g_seg=ggplot(pc_plot, aes(PC1,PC2)) +
  geom_point(size=3, color=NA) + ggtitle("Loadings for 48 DE genes") +
  scale_color_manual(values = c("red3", "mediumorchid3")) + coord_equal() + 
  geom_segment(data=datapc2, aes(x=0, y=0, xend=var1, yend=var2), 
               arrow=arrow(length=unit(0.2,"cm")), alpha=0.5) +
  geom_text_repel(data=datapc2, aes(x=var1, y=var2, label=varnames), size = 3)
g_seg

pdf("variable_genes_PCA_biplot_48DE.pdf", height = 8, width = 10)
grid.arrange(g_seg1+theme(legend.position = 'none'), g_seg, nrow = 1)
dev.off()



#
#
#
#
#
#
#
#
#
# Try a GOSeq analysis

library(goseq)
library(GenomicFeatures)
library("Mus.musculus")
seqlengths(Mus.musculus)

DEgenes <- read.table("Lab_v_wild_all.txt", T, '\t')

genes <- rep(0, nrow(DEgenes))
genes[DEgenes$log2FoldChange > 1.5 | DEgenes$log2FoldChange < -1.5 & DEgenes$padj < 0.05] <- 1
names(genes) <- DEgenes$Gene


pwf=nullp(genes,"mm10","geneSymbol")
head(pwf)

GO.wall=goseq(pwf,"mm10","geneSymbol")
GO.wall$padj <- p.adjust(GO.wall$over_represented_pvalue, method ="BH")
GO.sig <- subset(GO.wall, padj < 0.01)

write.table(GO.sig, "GOSeq/GO_enrichment.txt", row.names = F, sep='\t', quote = F)

DES <- names(genes[genes==1])
DES<- paste(DES, collapse=",")


GOTab <- data.frame(table(GO.wall$term))
