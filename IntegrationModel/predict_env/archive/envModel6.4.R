#### Predict Environment from multi-omic data types

### 4.8.19

setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/IntegrationModel/predict_env")

library(ggplot2)
library(caret)
library(plotROC)
library(randomForest)
library(matrixStats)


## read in metadata

meta <- readRDS("../ml_inputs/mouse_metadata.RData")

#
forester <- function(x_data, y_data, dataType, train_group, test_group) {
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
                             savePredictions = T, summaryFunction = twoClassSummary)
  ###500
  rf_test <- train(condition ~ ., 
                   data = data_run1[train_group,], tuneLength = 20, 
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
  
  g <- ggplot(huh, aes(m=wild, d=factor(obs, levels = c("lab", "wild")))) + 
    geom_roc(n.cuts=0) + 
    coord_equal() +
    style_roc()
  
  namer <- paste0(dataType, "_model_performance.pdf")
  pdf(namer, height = 5, width = 7)
  print(g + annotate("text", x=0.75, y=0.25, label=paste("AUC =", round((calc_auc(g))$AUC, 4))))
  dev.off()
  
  ##test
  #test results
  preds1 <- data.frame(X1=predict(rf_test, newdata = data_run1[test_group,]))
  #preds1$obs <- predict(lasso_test, newdata = data_test1)
  preds1$true <- data_run1[test_group,]$condition
  preds1$true <- gsub(0, "X0",preds1$true)
  preds1$true <- gsub(1, "X1",preds1$true)
  
  g <- ggplot(preds1, aes(m=X1, d=factor(true, levels = c("X0", "X1")))) + 
    geom_roc(n.cuts=0) + 
    coord_equal() +
    style_roc() +
    #ggtitle(paste0("Lasso_Test_", CV, "_", norm_names[f], "_", surv_num))
  glist[[counter]] <- g + annotate("text", x=0.75, y=0.25, label=paste("AUC =", round((calc_auc(g))$AUC, 4)))
  counter = counter+1
  
  #namer <- paste0(out_dir, "Lasso_Test_", CV, "_", surv_num, "_yrs_", norm_names[f], "_AUC.txt")
  #write.table(preds1, namer, sep='\t', row.names = F, quote = F)
  
  preds <- as.character(predict(rf_test, newdata = data_run1[test_group,]))
  cf <- confusionMatrix(data = factor(preds, levels=c("lab","wild")), 
                        reference = as.factor(data_run1[test_group,]$condition), 
                        mode = "prec_recall")
  
  c_mat <- c(truePos=cf$table[1,1], falsePos=cf$table[2,1], 
             falseNeg=cf$table[1,2], trueNeg=cf$table[2,2])
  
  #
  #
  #
  #
  ###log regress
  lr_test <- train(condition ~ ., 
                   data = data_run1[train_group,], tuneLength = 10, 
                   method = "regLogistic",
                   trControl = added_ctrl)
  
  namer <- paste0(dataType, "_rf_model.RData")
  saveRDS(lr_test, file=namer)
  
  # feature importance
  lr_imp <- varImp(lr_test)$importance
  
  
  bestModParam <- unlist(lr_test$bestTune)
  selectedIndices <- lr_test$pred == bestModParam
  huh <- lr_test$pred
  
  g <- ggplot(huh, aes(m=wild, d=factor(obs, levels = c("lab", "wild")))) + 
    geom_roc(n.cuts=0) + 
    coord_equal() +
    style_roc()
  
  namer <- paste0(dataType, "_model_performance.pdf")
  pdf(namer, height = 5, width = 7)
  print(g + annotate("text", x=0.75, y=0.25, label=paste("AUC =", round((calc_auc(g))$AUC, 4))))
  dev.off()
  
  ##test
  #test results
  preds1 <- data.frame(X1=predict(lr_test, newdata = data_run1[test_group,]))
  #preds1$obs <- predict(lasso_test, newdata = data_test1)
  preds1$true <- data_run1[test_group,]$condition
  preds1$true <- gsub(0, "X0",preds1$true)
  preds1$true <- gsub(1, "X1",preds1$true)
  
  g <- ggplot(preds1, aes(m=X1, d=factor(true, levels = c("lab", "wild")))) + 
    geom_roc(n.cuts=0) + 
    coord_equal() +
    style_roc() 
    #ggtitle(paste0("Lasso_Test_", CV, "_", norm_names[f], "_", surv_num))
    glist[[counter]] <- g + annotate("text", x=0.75, y=0.25, label=paste("AUC =", round((calc_auc(g))$AUC, 4)))
  counter = counter+1
  
  #namer <- paste0(out_dir, "Lasso_Test_", CV, "_", surv_num, "_yrs_", norm_names[f], "_AUC.txt")
  #write.table(preds1, namer, sep='\t', row.names = F, quote = F)
  
  preds <- as.character(predict(lr_test, newdata = data_run1[test_group,]))
  cf <- confusionMatrix(data = factor(preds, levels=c("lab","wild")), 
                        reference = as.factor(data_run1[test_group,]$condition), 
                        mode = "prec_recall")
  
  c_mat <- c(truePos=cf$table[1,1], falsePos=cf$table[2,1], 
             falseNeg=cf$table[1,2], trueNeg=cf$table[2,2])
}


#
#
#
### choose test and train groups
train_group <- as.character(sample(meta$mouse_id, ceiling(nrow(meta)*0.75)))
test_group <- meta$mouse_id[!meta$mouse_id %in% train_group]

### run the models
y_data <- meta$Environment
### rnaseq
x_data_genes <- readRDS("../ml_inputs/genes.RData")

# By DE
sig_genes <- as.character(unique(subset(
  read.table("../inputs/Lab_v_wild_all.txt", T, '\t'),
  abs(log2FoldChange) > 1.2 & padj < 0.05)$Gene))
x_data_genes <- x_data_genes[,sig_genes]

# By Variation
x_data_genes <- x_data_genes[,order(colVars(x_data_genes),decreasing=T)][,1:200]

# By Gene Module
x_data_genes <- readRDS("../ml_inputs/gene_modules.RData")

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

coul <- c("darkseagreen2", "lightcoral", "dodgerblue2", "mediumorchid3")
gs=c("genes", "cells", "cytokines", "otus")
color_coder <- c(Gene="lightcoral", ImmunePop="dodgerblue2", OTU="mediumorchid3", Cytokine="darkseagreen2")

feat_heat_filt <- feat_heat[order(feat_heat$Overall, decreasing=T),][1:50,]
feat_heat_filt$Features <- factor(feat_heat_filt$Features, levels = rev(feat_heat_filt$Features))
pdf("All_feature_imp_clean.pdf", height = 8, width = 6)
ggplot(feat_heat_filt, aes(Features, Overall, fill=DataType)) +
  geom_col() + coord_flip() + theme_bw() +
  ylab("MeanDecreaseGini") + ggtitle("Top 20 Model Features") +
  scale_fill_manual(values=color_coder)+
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

mat[mat<0.6 & mat > -0.6]=0
good_cols <- which(colSums(mat)!=1)

mat=mat[good_cols,good_cols]
idx=diag(ncol(mat))
mat[idx==1]<-0


library(ComplexHeatmap)
library(circlize)
Heatmap(mat,
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))

library(igraph)
#library(qgraph)
#network=graph_from_adjacency_matrix( mat, weighted=T, mode="undirected", diag=F)

#make facs names shorter
mln_cols <- grep("CD44", colnames(mat))
colnames(mat)[mln_cols] <- gsub("T_cells_", "T_cells\n", colnames(mat)[mln_cols])

mln_cols <- grep("MHCII", colnames(mat))
colnames(mat)[mln_cols] <- gsub("CD11c_", "CD11c\n", colnames(mat)[mln_cols])

#


##

network=graph_from_adjacency_matrix(mat, weighted=T, mode="undirected", diag=T)
icoord = layout_with_fr(network)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

E(network)$color[E(network)$weight > 0] <- add.alpha('red', alpha=0.4)
E(network)$color[E(network)$weight < 0] <- add.alpha('blue', alpha=0.4)

top_feats_meta$DataType <- factor(top_feats_meta$DataType, levels = c("Gene", "ImmunePop", 
                                                                      "OTU", "Cytokine"))
rownames(top_feats_meta) <- top_feats_meta$Features

my_color=color_coder[as.numeric(as.factor(top_feats_meta$DataType[as.numeric(good_cols)]))]
# plot

pdf("env_network6.4.pdf", height = 50, width = 50)
par(bg="white", mar=c(0,0,0,0))
set.seed(4)
plot(network, layout=icoord, 
     vertex.size=10,
     vertex.shape='circle',
     vertex.color=my_color, 
     vertex.label.cex=1.2,
     vertex.label.family="Helvetica",
     vertex.label.color="black",
     vertex.frame.color="transparent"
)
dev.off()


curr_coord = icoord
new_coord = data.frame(curr_coord, top_feats_meta[good_cols,])
colnames(new_coord) <- c("X", "Y", 'Feature', "Score", "DataType")

library(ggplot2)
ggplot(new_coord, aes(X,Y,color=DataType)) + 
  geom_point()

#gs = as.character(unique(new_coord$DataType))
gs=c("Gene", "ImmunePop", "OTU", "Cytokine")
better_coords <- data.frame(X=NA,Y=NA,Feature=NA,Score=NA,DataType=NA,X2=NA,Y2=NA)
low_bounds_y=c(250,800,250,0)
high_bounds_y=c(700,950,700,200)

low_bounds_x=c(0,300,700,300)
high_bounds_x=c(300,700,1000,700)

for(i in 1:length(gs)){
  curr <- subset(new_coord, DataType == gs[i])
  
  x2 = sample(low_bounds_x[i]:high_bounds_x[i], nrow(curr))
  y2 = sample(low_bounds_y[i]:high_bounds_y[i], nrow(curr), replace = T)
  
  curr$X2=x2
  curr$Y2=y2
  
  better_coords = rbind(better_coords, curr)
}
better_coords=better_coords[-1,]

ggplot(better_coords, aes(X2,Y2,color=DataType)) + 
  geom_point()
better_coords=better_coords[rownames(new_coord),]

pdf("env_network6.4.pdf", height = 50, width = 50)
par(bg="white", mar=c(0,0,0,0))
set.seed(4)
plot(network, layout=as.matrix(better_coords[,6:7]),
     vertex.size=10,
     vertex.color=my_color, 
     vertex.label.cex=1.2,
     vertex.label.family="Helvetica",
     vertex.label.color="black",
     vertex.frame.color="transparent"
)
dev.off()

tkplot(network, layout=as.matrix(better_coords[,6:7]),
       vertex.size=5,
       vertex.color=my_color, 
       vertex.label.cex=1.2,
       vertex.label.family="Helvetica",
       vertex.label.color="black",
       vertex.frame.color="grey90"
)


tk_coords = tk_coords(1)

pdf("env_network6.4.pdf", height = 40, width = 60)
par(bg="white", mar=c(0,0,0,0))
set.seed(4)
plot(network, #layout=as.matrix(better_coords[,5:6]),
     layout=tk_coords,
     vertex.size=10,
     vertex.color=my_color, 
     vertex.label.cex=1.2,
     vertex.label.family="Helvetica",
     vertex.label.color="black",
     vertex.frame.color="transparent"
)

legend(x=-1.5, y=-0.6, 
       legend=paste( levels(as.factor(top_feats_meta$DataType[good_cols]))), 
       col = coul[c(2,3,4,1)] , bty = "n", pch=20 , pt.cex = 14, cex = 9, text.col="black" , horiz = F)
dev.off()


#
#
#
#
#
#

### values
x_data_top <- data.frame(x_data_top)
x_data_top$Environment <- y_data

library(ggsignif)

ggplot(x_data_top, aes(Environment, log2(X9.2.2+1), color = Environment)) +
  geom_boxplot(alpha=0.6) + geom_jitter(width = 0.2) +
  scale_color_manual(values=c("mediumorchid", "red3")) +
  theme_bw() + geom_signif(test = "wilcox.test", 
                           comparisons = combn(levels(x_data_top$Environment),2, simplify = F),
                           color='black', map_signif_level = T)

ggplot(x_data_top, aes(Environment, log2(MLN_CD4_T_cells_CD25.+1), color = Environment)) +
  geom_boxplot(alpha=0.6, outlier.shape = NA) + geom_jitter(width = 0.2) +
  scale_color_manual(values=c("mediumorchid", "red3")) +
  theme_bw() + geom_signif(test = "wilcox.test", 
                           comparisons = combn(levels(x_data_top$Environment),2, simplify = F),
                           color='black', map_signif_level = T)

ggplot(x_data_top, aes(Environment, log2(IFN.y_CD3.CD28+1), color = Environment)) +
  geom_boxplot(alpha=0.6, outlier.shape = NA) + geom_jitter(width = 0.2) +
  scale_color_manual(values=c("mediumorchid", "red3")) +
  theme_bw() + geom_signif(test = "wilcox.test", 
                           comparisons = combn(levels(x_data_top$Environment),2, simplify = F),
                           color='black', map_signif_level = T)

x_data_sub <- x_data_top[c(1,4,23,51)]

library(gridExtra)
glist=list()
for(j in 1:3){
  curr <- x_data_sub[,c(j,4)]
  gg=ggplot(curr, aes(Environment, log2(curr[,1]+1), color = Environment)) +
    ggtitle(gsub("X","",colnames(curr)[1]))+
    ylab("")+
    geom_boxplot(alpha=0.6, outlier.shape = NA) + geom_jitter(width = 0.2) +
    scale_color_manual(values=c("mediumorchid", "red3")) +
    theme_bw() + geom_signif(test = "wilcox.test", 
                             comparisons = combn(levels(curr$Environment),2, simplify = F),
                             color='black', map_signif_level = T)
  glist[[j]]=grid.grabExpr(print(gg))
  
}


library(cowplot)
pdf("top_boxes.pdf", height = 4, width = 10)
grid.arrange(grobs=glist, nrow=1)
dev.off()

