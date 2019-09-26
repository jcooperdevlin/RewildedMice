######

##### We built our models

##### of top features what are groups in our data?

setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/IntegrationModel")

library(ggplot2)
library(ComplexHeatmap)
library(matrixStats)
library(caret)
library(gridExtra)
library(ggrepel)
library(circlize)

## read in metadata

meta <- readRDS("ml_inputs/mouse_metadata.RData")

### run the models
y_data <- meta$Environment
### rnaseq
x_data_genes <- readRDS("ml_inputs/genes.RData")

# By DE
#sig_genes <- as.character(unique(subset(
#  read.table("inputs/Lab_v_wild_all.txt", T, '\t'),
#  abs(log2FoldChange) > 1.2 & padj < 0.05)$Gene))
#x_data_genes <- x_data_genes[,sig_genes]

# By Variation
x_data_genes <- x_data_genes[,order(colVars(x_data_genes),decreasing=T)][,1:200]
x_data_genes <- scale(x_data_genes)

##otus
x_data_otus <- readRDS("ml_inputs/otus.RData")
x_data_otus <- scale(x_data_otus)

##facs
x_data_facs <- readRDS("ml_inputs/facs.RData")
x_data_facs <- scale(x_data_facs)

##cytokines
x_data_cyt <- readRDS("ml_inputs/cytokines.RData")
x_data_cyt <- scale(x_data_cyt)

#### all
x_data_all <- cbind(x_data_genes, x_data_otus, x_data_facs, x_data_cyt)


#### load models env
all_model <- readRDS("predict_env/All_rf_model.RData")
rf_imp <- data.frame(Overall = varImp(all_model)$importance$Overall)

good_gene_tab <- read.table("predict_env/All_feature_vector.txt", T, '\t')
x_data <- x_data_all
for(b in 1:ncol(x_data)){
  if(grepl(';',colnames(x_data)[b])){
    lookerupper <- paste0("^",colnames(x_data)[b], "$")
    if(grepl("\\[", lookerupper)){
      lookerupper1 <- strsplit(lookerupper, "\\[")
      lookerupper2 <- strsplit(lookerupper, "\\]")
      lookerupper <- paste0(lookerupper1[[1]][1], "\\[.*\\]",
                            lookerupper2[[1]][2])
    }
    lookupnum <- grep(lookerupper, good_gene_tab$Feature, fixed=F)
    colnames(x_data)[b] <- paste0("OTU_", lookupnum)
  }
}
feat_heat <- data.frame(Features = colnames(x_data))
feat_heat <- cbind(feat_heat, rf_imp)
feat_heat$DataType = c(rep("Gene", ncol(x_data_genes)),
                       rep("OTU", ncol(x_data_otus)),
                       rep("ImmunePop", ncol(x_data_facs)),
                       rep("Cytokine", ncol(x_data_cyt)))

feat_heat_env <- feat_heat[order(feat_heat$Overall, decreasing=T),]


###
#
#
#### load models env
all_model <- readRDS("predict_genotype//All_rf_model.RData")
rf_imp <- data.frame(Overall = varImp(all_model)$importance$Overall)

good_gene_tab <- read.table("predict_genotype//All_feature_vector.txt", T, '\t')
x_data <- x_data_all
for(b in 1:ncol(x_data)){
  if(grepl(';',colnames(x_data)[b])){
    lookerupper <- paste0("^",colnames(x_data)[b], "$")
    if(grepl("\\[", lookerupper)){
      lookerupper1 <- strsplit(lookerupper, "\\[")
      lookerupper2 <- strsplit(lookerupper, "\\]")
      lookerupper <- paste0(lookerupper1[[1]][1], "\\[.*\\]",
                            lookerupper2[[1]][2])
    }
    lookupnum <- grep(lookerupper, good_gene_tab$Feature, fixed=F)
    colnames(x_data)[b] <- paste0("OTU_", lookupnum)
  }
}
feat_heat <- data.frame(Features = colnames(x_data))
feat_heat <- cbind(feat_heat, rf_imp)
feat_heat$DataType = c(rep("Gene", ncol(x_data_genes)),
                       rep("OTU", ncol(x_data_otus)),
                       rep("ImmunePop", ncol(x_data_facs)),
                       rep("Cytokine", ncol(x_data_cyt)))

feat_heat_geno <- feat_heat[order(feat_heat$Overall, decreasing=T),]

#
#
#
#
#
## top 25 features from each model ~50 total and see if there are groups between env and genotype

combo_feats <- unique(c(as.character(feat_heat_env$Features[1:25]), as.character(feat_heat_geno$Features[1:25])))
combo_feats_df <- subset(rbind(feat_heat_env[1:25,], feat_heat_geno[1:25,]), Features %in% combo_feats)
combo_feats_df <- combo_feats_df[!duplicated(combo_feats_df$Features),]
rownames(combo_feats_df) <- as.character(combo_feats_df$Features)
combo_feats_df <- combo_feats_df[combo_feats,]

all_top <- x_data_all[,combo_feats]

pc <- prcomp(all_top, scale=F)
summary(pc)
expl_var <- pc$sdev^2/sum(pc$sdev^2)*100


plot.data = data.frame(meta, pc$x)


g0 <- ggplot(plot.data, aes(PC1, PC2, color = Environment, shape = Genotype)) +
  geom_point(size = 3) + 
  ylab(paste0("PC2 ", round(expl_var[2],2), "% expl. variation")) +
  xlab(paste0("PC1 ", round(expl_var[1],2), "% expl. variation")) +
  scale_color_manual(values = c("mediumpurple1", "red3")) +
  theme_bw()

grid.arrange(g0, g0 + facet_wrap(~ Genotype, ncol = 2), ncol = 1)


PC = pc
rownames(PC$x) <- 1:nrow(PC$x)
rownames(PC$rotation) <- colnames(all_top)
data <- data.frame(obsnames=row.names(PC$x), PC$x[,1:2])
datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation[,1:2])
datapc$var1 <- scales::rescale(datapc$PC1, c(min(data$PC1),max(data$PC1)))
datapc$var2 <- scales::rescale(datapc$PC2, c(min(data$PC2),max(data$PC2)))

datapc$mult <- abs(datapc$PC1*datapc$PC2)
datapc <- datapc[order(datapc$mult, decreasing = T),]
datapc2 = datapc[1:20,]

pc_seg=ggplot(plot.data, aes(PC1,PC2)) +
  geom_point(size=3, alpha=0) + #ggtitle("Loadings for Cytokine+Stim") +
  scale_color_manual(values = c("mediumorchid3", "red3")) + coord_equal() + 
  geom_segment(data=datapc2, aes(x=0, y=0, xend=var1, yend=var2), 
               arrow=arrow(length=unit(0.2,"cm")), alpha=0.5) +
  geom_text_repel(data=datapc2, aes(x=var1, y=var2, label=varnames), size = 3) +
  #ylim(c(-10,20))+xlim(c(-25,10))+
  theme_bw()


#
#
#

#
#
### make a heatmap

hh <- rowAnnotation(df = data.frame(Environment=meta$Environment, Genetics=meta$Genotype),
                    col = list(Environment = c(lab = "mediumorchid3",
                                               wild = "red3"),
                               Genetics = c(#AtgW = "darkorange1",
                                            AtgE = "dodgerblue2",
                                            AtgH = "navyblue",
                                            B6 = "darkorange1",
                                            NOD2 = "mediumspringgreen")))

ha <- columnAnnotation(df = data.frame(DataType=combo_feats_df$DataType),
                       col=list(DataType = c(ImmunePop="purple4",
                                             Gene="forestgreen",
                                             Cytokine="lightcoral",
                                             OTU="grey50")))
labs <- colnames(all_top)
ll <- columnAnnotation(labels = anno_text(labs, gp=gpar(fontsize=5),
                                          which = "column", rot=60, 
                                          just = 'right', offset=unit(3,"cm")), 
                       height = unit(1,"cm"))


pdf("top_features_heat.pdf", height = 5, width = 11)
Heatmap(all_top,
  show_row_names = F, show_column_names = F,
  top_annotation = ha,
  bottom_annotation = ll,
  bottom_annotation_height = unit(3,"cm"),
  heatmap_legend_param = list(title = "Feature\nValue"),
  cluster_rows = T, cluster_columns = T, #row_names_side = 'left',
  #row_names_gp = gpar(fontsize=8),
  #row_names_max_width = unit(10,'cm'),
  col = colorRamp2(c(-2,0,2), c("blue", "white", "red"))) + hh
dev.off()



