### Integration Model Test1
### 3.27

#setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/IntegrationModel")

library(ggplot2)
#library(MultivariateRandomForest)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(psych)
library(gridExtra)

source("glmnet_full_function.R")

## train data X
rnaseq_df <- read.table("normalizedCounts.txt", T, '\t')
rnaseq_var <- rnaseq_df[order(rowVars(data.matrix(rnaseq_df[,-1])),decreasing=T),]#[1:5000,]
trainX <- t(rnaseq_var[,-1])
colnames(trainX) <- rnaseq_var$GeneSymbol
rownames(trainX) <- gsub("X", "", rownames(trainX))
trainXs <- log2(trainX+1)

### train data y
cyt_df <- read.table("MLN_stimulation_flat.txt", T, '\t')
rownames(cyt_df) <- cyt_df$mouse_id
cyt_df <- log2(cyt_df[,-1]+1)
#cyt_df <- scale(cyt_df[,-1])
trainY <- data.matrix(cyt_df)
rownames(trainY) <- gsub("-", "_", rownames(trainY))

overlaps <- intersect(rownames(trainY), rownames(trainXs))
overlaps <- sample(overlaps, length(overlaps))

trainXs <- trainXs[overlaps,]
trainY <- trainY[overlaps,]

cytokines <- gsub("_Bac.*","",colnames(trainY))
cytokines <- gsub("_Cl.*","",cytokines)
cytokines <- gsub("_Ca.*","",cytokines)
cytokines <- gsub("_CD.*","",cytokines)
cytokines <- gsub("_Sta.*","",cytokines)
cytokines <- gsub("_P.*","",cytokines)

cytokines <- unique(cytokines)
for(i in 1:length(cytokines)){
  namer <- paste0("Gene_model/", cytokines[i])
  cyt_select <- grep(cytokines[i],colnames(trainY))
  run_glm(trainXs, trainY[,cyt_select], namer, type="Genes", p_threshold = 0.01)
}

#run_glm(trainXs, trainY, "Gene_model", type="Genes", p_threshold = 0.01)