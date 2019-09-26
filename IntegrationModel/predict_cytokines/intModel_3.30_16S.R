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



###### 16S Reads pvalue = 0.05

## train data X
otu_tab <- read.table("otu_table_counts.txt", T, '\t')
#otu_tab <- aggregate(otu_tab, by=list("taxonomy"), sum)
otu_tab <- aggregate(.~taxonomy, data=otu_tab, sum)
rownames(otu_tab) <- otu_tab$taxonomy
otu_tab <- otu_tab[,-1]
otu_tab <- t(otu_tab)
otu_tab <- otu_tab[,-1]

rownames(otu_tab) <- gsub("X", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("_A.*", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("_N.*", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("_B.*", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("_", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("\\.", "_", rownames(otu_tab))

trainXs <- log2(otu_tab+1)


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
  namer <- paste0("OTUs_counts_pval0.05/", cytokines[i])
  cyt_select <- grep(cytokines[i],colnames(trainY))
  run_glm(trainXs, trainY[,cyt_select], namer, type="OTU", p_threshold = 0.5)
}
#run_glm(trainXs, trainY, "OTUs_counts_pval0.05", type="OTUs", p_threshold = 0.05)




###### ###### 16S Reads pvalue = 0.1

## train data X
otu_tab <- read.table("otu_table_counts.txt", T, '\t')
otu_tab <- aggregate(.~taxonomy, data=otu_tab, sum)
rownames(otu_tab) <- otu_tab$taxonomy
otu_tab <- otu_tab[,-1]
otu_tab <- t(otu_tab)
otu_tab <- otu_tab[,-1]

rownames(otu_tab) <- gsub("X", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("_A.*", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("_N.*", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("_B.*", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("_", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("\\.", "_", rownames(otu_tab))

trainXs <- otu_tab


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
  namer <- paste0("OTUs_counts_pval0.5/", cytokines[i])
  cyt_select <- grep(cytokines[i],colnames(trainY))
  run_glm(trainXs, trainY[,cyt_select], namer, type="OTU", p_threshold = 0.5)
}

#run_glm(trainXs, trainY, "OTUs_counts_pval0.5", type="OTUs", p_threshold = 0.5)
#
#
#

#

###### Now train models with 16S relab pvalue=0.05

## train data X
otu_tab <- read.table("otu_table_relab.txt", T, '\t', skip=1)
rownames(otu_tab) <- otu_tab$ID
otu_tab <- otu_tab[,-1]
otu_tab <- t(otu_tab)
otu_tab <- otu_tab[,-1]

rownames(otu_tab) <- gsub("X", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("_A.*", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("_N.*", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("_B.*", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("_", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("\\.", "_", rownames(otu_tab))

trainXs <- otu_tab


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
  namer <- paste0("OTUs_relab_pval0.05/", cytokines[i])
  cyt_select <- grep(cytokines[i],colnames(trainY))
  run_glm(trainXs, trainY[,cyt_select], namer, type="OTU", p_threshold = 0.05)
}

#run_glm(trainXs, trainY, "OTUs_relab_pval0.05", type="OTUs", p_threshold = 0.05)



#
#

###### Now train models with 16S relab pvalue=0.05

## train data X
otu_tab <- read.table("otu_table_relab.txt", T, '\t', skip=1)
rownames(otu_tab) <- otu_tab$ID
otu_tab <- otu_tab[,-1]
otu_tab <- t(otu_tab)
otu_tab <- otu_tab[,-1]

rownames(otu_tab) <- gsub("X", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("_A.*", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("_N.*", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("_B.*", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("_", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("\\.", "_", rownames(otu_tab))

trainXs <- otu_tab


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
  namer <- paste0("OTUs_relab_pval0.5/", cytokines[i])
  cyt_select <- grep(cytokines[i],colnames(trainY))
  run_glm(trainXs, trainY[,cyt_select], namer, type="OTU", p_threshold = 0.5)
}

#run_glm(trainXs, trainY, "OTUs_relab_pval0.5", type="OTUs", p_threshold = 0.5)