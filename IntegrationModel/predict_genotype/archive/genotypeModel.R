#### Predict Environment from multi-omic data types

### 4.8.19

setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/IntegrationModel/predict_genotype")

library(ggplot2)
library(caret)
library(plotROC)
library(pROC)
library(randomForest)
library(MLmetrics)
library(matrixStats)


## read in metadata

meta <- readRDS("../ml_inputs/mouse_metadata.RData")

#
forester <- function(x_data, y_data, dataType) {
  good_gene_tab <- data.frame(num=1:length(colnames(x_data)), Feature=colnames(x_data))
  namer <- paste0(dataType, "_feature_vector.txt")
  write.table(good_gene_tab, namer, row.names=F, quote=F, sep='\t')
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
  
  data_run1 <- data.frame(x_data, condition=y_data)
  feat_heat <- data.frame(Features = colnames(x_data))
  
  added_ctrl <- trainControl(method = "repeatedcv", number = 5,
                             repeats = 10, verboseIter = T, 
                             classProbs = T,
                             savePredictions = T, summaryFunction = multiClassSummary)
  ###500
  rf_test <- train(condition ~ ., 
                   data = data_run1, tuneLength = 20, 
                   method = "rf",
                   trControl = added_ctrl)
  
  namer <- paste0(dataType, "_rf_model.RData")
  saveRDS(rf_test, file=namer)
  
  # feature importance
  rf_imp <- data.frame(Overall = varImp(rf_test)$importance$Overall)
  
  if(nrow(good_gene_tab)<50){
    nvar=nrow(good_gene_tab)
  } else {nvar=50}
  namer <- paste0(dataType, "_feature_imp.pdf")
  pdf(namer, height = 10, width = 7)
  varImpPlot(rf_test$finalModel, n.var=nvar)
  dev.off()
  
  feat_heat <- cbind(feat_heat, rf_imp)
  
  bestModParam <- unlist(rf_test$bestTune)
  selectedIndices <- rf_test$pred$mtry == bestModParam
  huh <- rf_test$pred[selectedIndices, ]
  
  #
  #
  
  #
  huh$new_obs <- as.numeric(huh$obs)
  huh$new_pred <- as.numeric(huh$pred)
  
  cc = multiclass.roc(huh$new_pred, huh$new_obs)$auc
  namer <- paste0(dataType, "_model_performance.txt")
  write(paste0("Multi-class AUC= ", round(cc,4)), namer)
}

### run the models
y_data <- meta$Genotype
### rnaseq
x_data_genes <- readRDS("../ml_inputs/genes.RData")
# By DE
sig_genes <- as.character(unique(subset(
  read.table("../inputs/Lab_v_wild_all.txt", T, '\t'),
  abs(log2FoldChange) > 1.2 & padj < 0.05)$Gene))
x_data_genes <- x_data_genes[,sig_genes]

# By Variation
x_data_genes <- x_data_genes[,order(colVars(x_data_genes),decreasing=T)][,1:200]


forester(x_data = x_data_genes, y_data = y_data, dataType = "Gene")

##otus
x_data_otus <- readRDS("../ml_inputs/otus.RData")
forester(x_data = x_data_otus, y_data = y_data, dataType = "OTU")

##facs
x_data_facs <- readRDS("../ml_inputs/facs.RData")
forester(x_data = x_data_facs, y_data = y_data, dataType = "FACS")

##cytokines
x_data_cyt <- readRDS("../ml_inputs/cytokines.RData")
forester(x_data = x_data_cyt, y_data = y_data, dataType = "Cytokine")


##
#
#

#### all
x_data_all <- cbind(x_data_genes, x_data_otus, x_data_facs, x_data_cyt)
forester(x_data = x_data_all, y_data = y_data, dataType = "All")

#
#

all_model <- readRDS("All_rf_model.RData")
rf_imp <- data.frame(Overall = varImp(all_model)$importance$Overall)

good_gene_tab <- read.table("All_feature_vector.txt", T, '\t')
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

feat_heat_filt <- feat_heat[order(feat_heat$Overall, decreasing=T),][1:20,]
feat_heat_filt$Features <- factor(feat_heat_filt$Features, levels = rev(feat_heat_filt$Features))
pdf("All_feature_imp_clean.pdf", height = 4, width = 6)
ggplot(feat_heat_filt, aes(Features, Overall, fill=DataType)) +
  geom_col() + coord_flip() + theme_bw() +
  ylab("MeanDecreaseGini") + ggtitle("Top 20 Model Features") +
  scale_fill_manual(values=c(ImmunePop="dodgerblue2", 
                             Gene="forestgreen", 
                             Cytokine="darkorange2"))+
  theme(axis.text.y = element_text(size=12, color="black"),
        axis.text.x = element_blank(),
        axis.title = element_text(size=13, color="black"),
        legend.text = element_text(size=12, color="black"),
        legend.title = element_text(size=13, color="black"),
        plot.title = element_text(size=14, face='bold'))
dev.off()


#
#
#
#
## network?
top_feats <- as.character(feat_heat[order(feat_heat$Overall, decreasing=T),][1:50,]$Features)
top_feats_meta <- feat_heat[order(feat_heat$Overall, decreasing=T),][1:50,]

x_data_top <- x_data[,top_feats]
mat=cor(x_data_top)
mat[mat<0.65 & mat > -0.65]=0
good_cols <- which(colSums(mat)>1)

mat=mat[good_cols,good_cols]

#Heatmap(mat,
#        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))

library(igraph)
#library(qgraph)
network=graph_from_adjacency_matrix( mat, weighted=T, mode="undirected", diag=F)

coul <- c("darkorange2", "forestgreen", "dodgerblue2", "purple4")
my_color=coul[as.numeric(as.factor(top_feats_meta$DataType[good_cols]))]

#icoord = layout_with_kk(network)
#icoord = layout_with_drl(network)
#icoord = layout_with_dh(network)
#e <- get.edgelist(network)
#l <- qgraph.layout.fruchtermanreingold(e, vcount=vcount(network))

# plot
E(network)$color[E(network)$weight > 0] <- 'red'
E(network)$color[E(network)$weight < 0] <- 'blue'

pdf("genotype_top_feats_network.pdf", height = 10, width = 15)
par(bg="white", mar=c(0,0,0,0))
set.seed(4)
tkplot(network, #layout=icoord, 
     vertex.size=12,
     vertex.color=my_color, 
     vertex.label.cex=1.3,
     vertex.label.family="Helvetica",
     vertex.label.color="black"#,
     #vertex.frame.color="transparent"
)
# title and legend
legend(x=-1, y=-0.6, 
       legend=paste( levels(as.factor(top_feats_meta$DataType))), 
       col = coul , bty = "n", pch=20 , pt.cex = 2, cex = 1, text.col="black" , horiz = F)
dev.off()

