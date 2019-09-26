### supplementary table prep 8.5

setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/")


### ST 2 metadata FACs and 

meta <- read.table("int/data/metadata/mice_metadata.11.19_mouse_id.txt", T, '\t')


mln_cyt <- read.table("int/data/MLN_stimulations/MLN_stimulation_flat.txt", T, "\t")
plasma_cyt <- read.table("int/data/Serum_stimulations/Serum_names_final.txt", T, "t")
plasma_cyt <- read.table("int/data/Serum_stimulations/Serum_names_final.txt", T, "\t")

meta_plasma <- merge(meta, plasma_cyt, by.x = "mouse_id", by.y = "sample", all.x=T)

meta_plasma_cyt <- merge(meta_plasma, mln_cyt, by.x = "mouse_id", by.y = "mouse_id", all.x=T)

colnames(meta_plasma_cyt)

#fix bugs
colnames(meta_plasma_cyt) <- gsub("CD3\\.CD28", "aCD3-CD28", colnames(meta_plasma_cyt))
colnames(meta_plasma_cyt) <- gsub("\\.", "-", colnames(meta_plasma_cyt))
colnames(meta_plasma_cyt) <- gsub("ClostridiumP", "C. perfringens", colnames(meta_plasma_cyt))
colnames(meta_plasma_cyt) <- gsub("PseudomonasA", "P. aeruginosa", colnames(meta_plasma_cyt))
colnames(meta_plasma_cyt) <- gsub("BacillusS", "B. subtilis", colnames(meta_plasma_cyt))
colnames(meta_plasma_cyt) <- gsub("StaphA", "S. aureus", colnames(meta_plasma_cyt))
colnames(meta_plasma_cyt) <- gsub("BacteroidesV", "B. vulgatus", colnames(meta_plasma_cyt))
colnames(meta_plasma_cyt) <- gsub("CandidaA", "C. albicans", colnames(meta_plasma_cyt))
colnames(meta_plasma_cyt) <- gsub("_", " ", colnames(meta_plasma_cyt))



### FACS
name_change <- read.table("int/data/FACS_data/lymph_name_change.txt", F, '\t')
blood_lymph_facs <- read.table("int/data/FACS_data/BLOOD_lymph_FACS_metadata_names.txt", T, '\t')
mln_lymph_facs <- read.table("int/data/FACS_data/MLN_lymph_FACS_metadata_names.txt", T, '\t')
mln_myeloid_facs <- read.table("int/data/FACS_data/MLN_myeloid_FACS_metadata_names.txt", T, '\t')

rownames(blood_lymph_facs) <- blood_lymph_facs$mouse_id
colnames(blood_lymph_facs)[13:ncol(blood_lymph_facs)] <- as.character(name_change$V2)
colnames(blood_lymph_facs)[13:ncol(blood_lymph_facs)] <- paste0("Blood_", colnames(blood_lymph_facs)[13:ncol(blood_lymph_facs)])

rownames(mln_lymph_facs) <- mln_lymph_facs$mouse_id
colnames(mln_lymph_facs)[13:ncol(blood_lymph_facs)] <- as.character(name_change$V2)
colnames(mln_lymph_facs)[13:ncol(mln_lymph_facs)] <- paste0("MLN_", colnames(mln_lymph_facs)[13:ncol(mln_lymph_facs)])

rownames(mln_myeloid_facs) <- mln_myeloid_facs$mouse_id
colnames(mln_myeloid_facs)[13:ncol(mln_myeloid_facs)] <- paste0("MLN_", colnames(mln_myeloid_facs)[13:ncol(mln_myeloid_facs)])

lymph <- merge(blood_lymph_facs[,c(1,13:ncol(blood_lymph_facs))],
               mln_lymph_facs[,c(1,13:ncol(mln_lymph_facs))],
               by = "mouse_id")
full_facs <- merge(lymph, mln_myeloid_facs[,c(1,13:ncol(mln_myeloid_facs))],
                   by = "mouse_id")
colnames(full_facs) <- gsub("_", " ", colnames(full_facs))

full_facs$`mouse id` <- factor(c("02-3", "02-4", "02-7", "02-8", "02-9", as.character(full_facs$`mouse id`[6:180])))

##

meta_plasma_cyt_facs <- merge(meta_plasma_cyt, full_facs, by.x="mouse id", by.y="mouse id", all.x=T)
rownames(meta_plasma_cyt_facs) <- as.character(meta_plasma_cyt_facs$`mouse id`)
work=as.character(read.table("int/data/metadata/workidid", F, '\t')$V1)
meta_plasma_cyt_facs <- meta_plasma_cyt_facs[work,]
write.table(meta_plasma_cyt_facs, "int/data/metadata/ST2.txt", row.names=F, quote=F, sep='\t')


##### ST3 supervised int model data table

setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/IntegrationModel/predict_env")

meta <- readRDS("../ml_inputs/mouse_metadata.RData")

### run the models
names <- meta$mouse_id
y_data1 <- as.character(meta$Environment)
y_data2 <- as.character(meta$Genotype)

### rnaseq
x_data_genes <- readRDS("../ml_inputs/genes.RData")
# genes
x_data_genes <- x_data_genes[,order(colVars(x_data_genes),decreasing=T)][,1:200]
##otus
x_data_otus <- readRDS("../ml_inputs/otus.RData")
##facs
x_data_facs <- readRDS("../ml_inputs/facs.RData")
##cytokines
x_data_cyt <- readRDS("../ml_inputs/cytokines.RData")


#### all
x_data_all <- cbind(x_data_genes, x_data_otus, x_data_facs, x_data_cyt)
colnames(x_data_all) <- gsub("CD3\\.CD28", "aCD3-CD28", colnames(x_data_all))
colnames(x_data_all) <- gsub("\\.", "-", colnames(x_data_all))
colnames(x_data_all) <- gsub("ClostridiumP", "C. perfringens", colnames(x_data_all))
colnames(x_data_all) <- gsub("PseudomonasA", "P. aeruginosa", colnames(x_data_all))
colnames(x_data_all) <- gsub("BacillusS", "B. subtilis", colnames(x_data_all))
colnames(x_data_all) <- gsub("StaphA", "S. aureus", colnames(x_data_all))
colnames(x_data_all) <- gsub("BacteroidesV", "B. vulgatus", colnames(x_data_all))
colnames(x_data_all) <- gsub("CandidaA", "C. albicans", colnames(x_data_all))
colnames(x_data_all)[437:576] <- gsub("_", " ", colnames(x_data_all)[437:576])

setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/")
write.table(cbind(names, y_data1, y_data2, x_data_all), "int/data/metadata/ST3.txt", row.names=F, quote=F, sep='\t')


#
#
#
#

##### ST4 supervised int model data table

setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/IntegrationModel/network/")

meta <- readRDS("../ml_inputs/mouse_metadata.RData")

### run the models
names <- meta$mouse_id
y_data1 <- as.character(meta$Environment)
y_data2 <- as.character(meta$Genotype)

genes = readRDS("../ml_inputs/gene_modules.RData")
cytokines = readRDS("../ml_inputs/cytokines.RData")
cells = readRDS("../ml_inputs/facs.RData")
#otus = readRDS("../ml_inputs/otus_relab.RData")
otus = readRDS("../ml_inputs/otus.RData")

#### all
x_data_all <- cbind(genes, otus, cells, cytokines)
colnames(x_data_all) <- gsub("CD3\\.CD28", "aCD3-CD28", colnames(x_data_all))
colnames(x_data_all)[447:586] <- gsub("\\.", "-", colnames(x_data_all))[447:586]
colnames(x_data_all) <- gsub("ClostridiumP", "C. perfringens", colnames(x_data_all))
colnames(x_data_all) <- gsub("PseudomonasA", "P. aeruginosa", colnames(x_data_all))
colnames(x_data_all) <- gsub("BacillusS", "B. subtilis", colnames(x_data_all))
colnames(x_data_all) <- gsub("StaphA", "S. aureus", colnames(x_data_all))
colnames(x_data_all) <- gsub("BacteroidesV", "B. vulgatus", colnames(x_data_all))
colnames(x_data_all) <- gsub("CandidaA", "C. albicans", colnames(x_data_all))
colnames(x_data_all)[447:586] <- gsub("_", " ", colnames(x_data_all)[447:586])

setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/")
write.table(cbind(names, y_data1, y_data2, x_data_all), "int/data/metadata/ST4.txt", row.names=F, quote=F, sep='\t')

