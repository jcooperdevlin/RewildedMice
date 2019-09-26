#### Predict Genotype from multi-omic data types

### 6.21.19

setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/IntegrationModel/predict_genotype/")

library(ggplot2)
library(caret)
library(plotROC)
library(randomForest)
library(matrixStats)
library(pROC)


## read in metadata

meta <- readRDS("../ml_inputs/mouse_metadata.RData")

#modeler(x_data = x_data_genes, y_data = y_data, dataType = "Gene", 
#        train_group, test_group, "regLogistic")

#
modeler <- function(x_data, y_data, dataType, train_group, test_group, model) {
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
  mod_test <- train(condition ~ ., 
                   data = data_run1[train_group,], tuneLength = 10, 
                   method = model,
                   trControl = added_ctrl,
                   MaxNWts =10000000)
  
  namer <- paste0(dataType, "_", model, "_model.RData")
  saveRDS(mod_test, file=namer)
  
  # feature importance
  imp_df <- varImp(mod_test)$importance
  feat_heat <- cbind(feat_heat, imp_df)
  namer <- paste0(dataType, "_", model, "_featImp.txt")
  write.table(feat_heat, namer, quote=F, sep='\t', row.names=F)
  
  
  predictions <- as.numeric(predict(mod_test, data_run1[test_group,]))
  roc.multi <- multiclass.roc(as.numeric(data_run1[test_group,]$condition), predictions)
  cc=auc(roc.multi)
  
  namer <- paste0(dataType, "_", model, "_performance.txt")
  write(paste0("Multi-class AUC= ", round(cc,4)), namer)
  
  preds <- as.character(predict(mod_test, newdata = data_run1[test_group,]))
  cf <- confusionMatrix(data = factor(preds, levels=c("AtgH","AtgE", "NOD2", "B6")), 
                        reference = as.factor(data_run1[test_group,]$condition), 
                        mode = "prec_recall")
  
  c_mat <- cf$table
  return(c_mat)

}
#
#
#
#
#
#

### choose test and train groups
train_group <- as.character(sample(meta$mouse_id, ceiling(nrow(meta)*0.75)))
test_group <- meta$mouse_id[!meta$mouse_id %in% train_group]
test_meta <- subset(meta, mouse_id %in% test_group)
table(test_meta$Genotype)

### run the models
y_data <- meta$Genotype
### rnaseq
x_data_genes <- readRDS("../ml_inputs/genes.RData")

# By DE
#sig_genes <- as.character(unique(subset(
#  read.table("../inputs/Lab_v_wild_all.txt", T, '\t'),
#  abs(log2FoldChange) > 1.2 & padj < 0.05)$Gene))
#x_data_genes <- x_data_genes[,sig_genes]

# By Variation
x_data_genes <- x_data_genes[,order(colVars(x_data_genes),decreasing=T)][,1:200]

# By Gene Module
#x_data_genes <- readRDS("../ml_inputs/gene_modules.RData")

modeler(x_data = x_data_genes, y_data = y_data, dataType = "Gene", 
        train_group, test_group, "rf")
modeler(x_data = x_data_genes, y_data = y_data, dataType = "Gene", 
        train_group, test_group, "multinom")




##otus
x_data_otus <- readRDS("../ml_inputs/otus.RData")
modeler(x_data = x_data_otus, y_data = y_data, dataType = "OTU", 
        train_group, test_group, "rf")
modeler(x_data = x_data_otus, y_data = y_data, dataType = "OTU", 
        train_group, test_group, "multinom")

##facs
x_data_facs <- readRDS("../ml_inputs/facs.RData")
modeler(x_data = x_data_facs, y_data = y_data, dataType = "FACS", 
        train_group, test_group, "rf")
modeler(x_data = x_data_facs, y_data = y_data, dataType = "FACS", 
        train_group, test_group, "multinom")

##cytokines
x_data_cyt <- readRDS("../ml_inputs/cytokines.RData")
modeler(x_data = x_data_cyt, y_data = y_data, dataType = "Cytokine", 
        train_group, test_group, "rf")
modeler(x_data = x_data_cyt, y_data = y_data, dataType = "Cytokine", 
        train_group, test_group, "multinom")


##
#
#

#### all
x_data_all <- cbind(x_data_genes, x_data_otus, x_data_facs, x_data_cyt)
modeler(x_data = x_data_all, y_data = y_data, dataType = "All", 
        train_group, test_group, "rf")
modeler(x_data = x_data_all, y_data = y_data, dataType = "All", 
        train_group, test_group, "multinom")
