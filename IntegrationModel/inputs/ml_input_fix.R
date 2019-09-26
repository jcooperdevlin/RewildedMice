### clean up inputs

### 4.8.19
library(matrixStats)

setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/IntegrationModel/inputs")

#### Data types = MLN RNA-Seq, FACS, 16S, MLN Cytokines, Serum Cytokines

#
#
#
##### RNA-Seq

## train data RNASeq
rnaseq_df <- read.table("normalizedCounts.txt", T, '\t')
rnaseq_var <- rnaseq_df[order(rowVars(data.matrix(rnaseq_df[,-1])),decreasing=T),]#[1:5000,]
trainX <- t(rnaseq_var[,-1])
colnames(trainX) <- rnaseq_var$GeneSymbol
rownames(trainX) <- gsub("X", "", rownames(trainX))
#trainXs_genes <- log2(trainX+1)
trainXs_genes <- trainX

#or gene modules
gene_modules_df <- read.table("gene_modules_210.txt", T, '\t')
trainX <- t(gene_modules_df[,-1])
colnames(trainX) <- gene_modules_df$Module
rownames(trainX) <- gsub("X", "", rownames(trainX))
#trainXs_genes <- log2(trainX+1)
trainXs_modules <- trainX



## train data X FACS
name_change <- read.table("lymph_name_change.txt", F, '\t')

blood_lymph_facs <- read.table("BLOOD_lymph_FACS_metadata_names.txt", T, sep='\t')
rownames(blood_lymph_facs) <- blood_lymph_facs$mouse_id
colnames(blood_lymph_facs)[13:ncol(blood_lymph_facs)] <- as.character(name_change$V2)
colnames(blood_lymph_facs)[13:ncol(blood_lymph_facs)] <- paste0("Blood_", colnames(blood_lymph_facs)[13:ncol(blood_lymph_facs)])

mln_lymph_facs <- read.table("MLN_lymph_FACS_metadata_names.txt", T, sep='\t')
rownames(mln_lymph_facs) <- mln_lymph_facs$mouse_id
colnames(mln_lymph_facs)[13:ncol(blood_lymph_facs)] <- as.character(name_change$V2)
colnames(mln_lymph_facs)[13:ncol(mln_lymph_facs)] <- paste0("MLN_", colnames(mln_lymph_facs)[13:ncol(mln_lymph_facs)])

mln_myeloid_facs <- read.table("MLN_myeloid_FACS_metadata_names.txt", T, sep='\t')
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

#trainXs_otus <- log2(otu_tab+1)
trainXs_otus <- otu_tab


##otu relab

## train data X OTU table
otu_tab <- read.table("otu_table_relab.txt", T, '\t')

otu_tab <- aggregate(.~ID, data=otu_tab, sum)
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

#trainXs_otus <- log2(otu_tab+1)
trainXs_otus_relab <- otu_tab



### train data y
cyt_df <- read.table("MLN_stimulation_flat.txt", T, '\t')
rownames(cyt_df) <- cyt_df$mouse_id
cyt_df <- cyt_df[,-1]
#cyt_df <- log2(cyt_df[,-1]+1)
trainXs_cyt <- data.matrix(cyt_df)
rownames(trainXs_cyt) <- gsub("-", "_", rownames(trainXs_cyt))

#
#
#
#ser_df <- read.table("Serum_names_final.txt", T, '\t')
#rownames(ser_df) <- ser_df$sample
#ser_df <- log2(ser_df[,-1]+1)
#ser_df <- ser_df[,-1]
#trainXs_ser <- data.matrix(ser_df)
#rownames(trainXs_ser) <- gsub("-", "_", rownames(trainXs_ser))
#
#so few overlaps serum cannot be added to model

all_overlaps <- Reduce(intersect, list(rownames(trainXs_genes),
                                       rownames(trainXs_modules),
                                       rownames(trainXs_facs),
                                       rownames(trainXs_otus),
                                       rownames(trainXs_otus_relab),
                                       rownames(trainXs_cyt)))


saveRDS(trainXs_genes[all_overlaps,], file="../ml_inputs/genes.RData")
saveRDS(trainXs_modules[all_overlaps,], file="../ml_inputs/gene_modules.RData")
saveRDS(trainXs_facs[all_overlaps,], file="../ml_inputs/facs.RData")
saveRDS(trainXs_otus[all_overlaps,], file="../ml_inputs/otus.RData")
saveRDS(trainXs_otus_relab[all_overlaps,], file="../ml_inputs/otus_relab.RData")
saveRDS(trainXs_cyt[all_overlaps,], file="../ml_inputs/cytokines.RData")


#
#

#
#
### fix metadata

meta <- read.table("mice_metadata.11.19_mouse_id.txt", T, '\t')
meta$mouse_id <- gsub("-", "_", meta$mouse_id)
meta <- subset(meta, mouse_id %in% all_overlaps)
rownames(meta) <- meta$mouse_id
meta<-meta[all_overlaps,]
saveRDS(meta, file= "../ml_inputs/mouse_metadata.RData")
