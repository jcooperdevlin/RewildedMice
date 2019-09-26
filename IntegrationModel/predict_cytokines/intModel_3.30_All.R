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

## train data X RNASeq
rnaseq_df <- read.table("normalizedCounts.txt", T, '\t')
rnaseq_var <- rnaseq_df[order(rowVars(data.matrix(rnaseq_df[,-1])),decreasing=T),]#[1:5000,]
trainX <- t(rnaseq_var[,-1])
colnames(trainX) <- rnaseq_var$GeneSymbol
rownames(trainX) <- gsub("X", "", rownames(trainX))
trainXs_genes <- log2(trainX+1)

## train data X FACS
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
trainXs_facs <- data.matrix(trainXs)

## train data X OTU table
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

trainXs_otus <- log2(otu_tab+1)

all_overlaps <- intersect(
  intersect(rownames(trainXs_genes), rownames(trainXs_facs)),
  rownames(trainXs_otus))

trainXs <- cbind(trainXs_genes[all_overlaps,],
                 trainXs_facs[all_overlaps,],
                 trainXs_otus[all_overlaps,])


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
  namer <- paste0("All_model/", cytokines[i])
  cyt_select <- grep(cytokines[i],colnames(trainY))
  run_glm(trainXs, trainY[,cyt_select], namer, type="All", p_threshold = 0.05)
}

#run_glm(trainXs, trainY, "All_model", type="All", p_threshold = 0.05)
