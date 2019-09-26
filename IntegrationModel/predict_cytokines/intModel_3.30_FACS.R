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


###### Now train models with FACS data 

## train data X
blood_lymph_facs <- read.table("BLOOD_lymph_FACS_metadata_names.txt", T, sep='\t')
remove <- which(duplicated(blood_lymph_facs$mouse_id))
blood_lymph_facs <- blood_lymph_facs[-remove,]
rownames(blood_lymph_facs) <- blood_lymph_facs$mouse_id
colnames(blood_lymph_facs)[13:ncol(blood_lymph_facs)] <- paste0("Blood_", colnames(blood_lymph_facs)[13:ncol(blood_lymph_facs)])

mln_lymph_facs <- read.table("MLN_lymph_FACS_metadata_names.txt", T, sep='\t')
remove <- which(duplicated(mln_lymph_facs$mouse_id))
mln_lymph_facs <- mln_lymph_facs[-remove,]
rownames(mln_lymph_facs) <- mln_lymph_facs$mouse_id
colnames(mln_lymph_facs)[13:ncol(mln_lymph_facs)] <- paste0("MLN_", colnames(mln_lymph_facs)[13:ncol(mln_lymph_facs)])

mln_myeloid_facs <- read.table("MLN_myeloid_FACS_metadata_names.txt", T, sep='\t')
remove <- which(duplicated(mln_myeloid_facs$mouse_id))
mln_myeloid_facs <- mln_myeloid_facs[-remove,]
rownames(mln_myeloid_facs) <- mln_myeloid_facs$mouse_id
colnames(mln_myeloid_facs)[13:ncol(mln_myeloid_facs)] <- paste0("MLN_", colnames(mln_myeloid_facs)[13:ncol(mln_myeloid_facs)])


lymph <- merge(blood_lymph_facs[,c(1,13:ncol(blood_lymph_facs))],
               mln_lymph_facs[,c(1,13:ncol(mln_lymph_facs))],
               by = "mouse_id")
full_facs <- merge(lymph, mln_myeloid_facs[,c(1,13:ncol(mln_myeloid_facs))],
                   by = "mouse_id")
remover <- unique(which(is.na(full_facs), arr.ind = T)[,1])
full_facs <- full_facs[-remover,]
trainXs <- full_facs[,-1]
rownames(trainXs) <- gsub("-", "_", full_facs$mouse_id)
trainXs <- data.matrix(trainXs)


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
  namer <- paste0("FACS_model/", cytokines[i])
  cyt_select <- grep(cytokines[i],colnames(trainY))
  run_glm(trainXs, trainY[,cyt_select], namer, type="FACS", p_threshold = 0.05)
}


