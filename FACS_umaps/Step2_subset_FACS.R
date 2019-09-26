#### With UMAP or TSNE coordinates plot distributions and subset by major markers 

### Blood
umap_blood_df <- read.table("inputs/umap_combo_Blood.txt", F, '\t')

meta <- read.table("mice_metadata.11.19_mouse_id.txt", T, '\t')
orig_df <- read.table("inputs/Blood_df.txt", F, '\t')
ids <- orig_df$V1
names <- colnames(read.table("name_change.csv", T, ","))
colnames(orig_df) <- c("id", names)

### add metadata

orig_df$id <- factor(orig_df$id, levels = unique(orig_df$id))
orig_df <- orig_df[order(orig_df$id),]

rownames(meta) <- meta$mouse_id
meta <- meta[levels(orig_df$id),]

uniq_ids <- unique(orig_df$id)
orig_df2 <- data.frame(orig_df[1,], Genotype = NA, Environment = NA,
                       Wedge_cage = NA, Gender = NA, Pregnant = NA, Diarrhea = NA, Flow.date=NA)
for (j in 1:length(uniq_ids)){
  curr <- subset(orig_df, id == uniq_ids[j])
  meta_curr <- subset(meta, mouse_id == as.character(uniq_ids[j]))
  
  curr$Genotype <- rep(meta_curr$Genotype, each = nrow(curr))
  curr$Environment <- rep(meta_curr$Environment, each =nrow(curr))
  curr$Wedge_cage <- rep(meta_curr$Wedge_cage, each =nrow(curr))
  curr$Gender <- rep(meta_curr$Gender, each =nrow(curr))
  curr$Pregnant <- rep(meta_curr$Pregnant, each =nrow(curr))
  curr$Diarrhea <- rep(meta_curr$Diarrhea, each =nrow(curr))
  curr$Flow.date <- rep(meta_curr$Flow.date, each =nrow(curr))
  
  orig_df2 <- rbind(orig_df2, curr)
}
orig_df <- orig_df2[-1,]

colnames(orig_df)

# blood
orig_df$umap_1 <- umap_blood_df$V1
orig_df$umap_2 <- umap_blood_df$V2

umap_xlims <- c(floor(min(orig_df$umap_1)), ceiling(max(orig_df$umap_1)))
umap_ylims <- c(floor(min(orig_df$umap_2)), ceiling(max(orig_df$umap_2)))

namer <- paste0("Distributions/Blood_distr.pdf")
pdf(namer, height=15,width=10)
par(mfrow = c(4, 4))
for(i in 2:17){
  
  ylims=c(0,0.01)
  plot(density(orig_df[,i]), 
       main=paste0(colnames(orig_df)[i]),ylim=ylims)
}
dev.off()

thresh <- read.table("Distributions/thresholds/Blood_major_thresholds.txt", T, "\t")

orig_df$CD3[orig_df$CD3<thresh$CD3[1]] <- 0
orig_df$CD19[orig_df$CD19<thresh$CD19[1]] <- 0
orig_df$CD4[orig_df$CD4<thresh$CD4[1]] <- 0
orig_df$CD8[orig_df$CD8<thresh$CD8[1]] <- 0

#### write new df of CD4+ CD8+ CD19+ and CD3+
og_df <- read.table("inputs/Blood_df.txt", F, '\t')

#CD3
write.table(og_df[which(orig_df$CD3>thresh$CD3[1]),], "inputs/CD3_Blood_df.txt", 
            sep='\t', row.names=F, quote=F, col.names = F)

#CD19
write.table(og_df[which(orig_df$CD19>thresh$CD19[1]),], "inputs/CD19_Blood_df.txt", 
            sep='\t', row.names=F, quote=F, col.names = F)

#CD4
write.table(og_df[which(orig_df$CD4>thresh$CD4[1]),], "inputs/CD4_Blood_df.txt", 
            sep='\t', row.names=F, quote=F, col.names = F)

#CD8
write.table(og_df[which(orig_df$CD8>thresh$CD8[1]),], "inputs/CD8_Blood_df.txt", 
            sep='\t', row.names=F, quote=F, col.names = F)


#
#
#
#
#
#
### MLN
umap_mln_df <- read.table("inputs/umap_combo_MLN.txt", F, '\t')

meta <- read.table("mice_metadata.11.19_mouse_id.txt", T, '\t')
orig_df <- read.table("inputs/MLN_df.txt", F, '\t')
ids <- orig_df$V1
names <- colnames(read.table("name_change.csv", T, ","))
colnames(orig_df) <- c("id", names)

### add metadata

orig_df$id <- factor(orig_df$id, levels = unique(orig_df$id))
orig_df <- orig_df[order(orig_df$id),]

rownames(meta) <- meta$mouse_id
meta <- meta[levels(orig_df$id),]

uniq_ids <- unique(orig_df$id)
orig_df2 <- data.frame(orig_df[1,], Genotype = NA, Environment = NA,
                       Wedge_cage = NA, Gender = NA, Pregnant = NA, Diarrhea = NA, Flow.date=NA)
for (j in 1:length(uniq_ids)){
  curr <- subset(orig_df, id == uniq_ids[j])
  meta_curr <- subset(meta, mouse_id == as.character(uniq_ids[j]))
  
  curr$Genotype <- rep(meta_curr$Genotype, each = nrow(curr))
  curr$Environment <- rep(meta_curr$Environment, each =nrow(curr))
  curr$Wedge_cage <- rep(meta_curr$Wedge_cage, each =nrow(curr))
  curr$Gender <- rep(meta_curr$Gender, each =nrow(curr))
  curr$Pregnant <- rep(meta_curr$Pregnant, each =nrow(curr))
  curr$Diarrhea <- rep(meta_curr$Diarrhea, each =nrow(curr))
  curr$Flow.date <- rep(meta_curr$Flow.date, each =nrow(curr))
  
  orig_df2 <- rbind(orig_df2, curr)
}
orig_df <- orig_df2[-1,]

colnames(orig_df)

# mln
orig_df$umap_1 <- umap_mln_df$V1
orig_df$umap_2 <- umap_mln_df$V2

umap_xlims <- c(floor(min(orig_df$umap_1)), ceiling(max(orig_df$umap_1)))
umap_ylims <- c(floor(min(orig_df$umap_2)), ceiling(max(orig_df$umap_2)))

namer <- paste0("Distributions/MLN_distr.pdf")
pdf(namer, height=15,width=10)
par(mfrow = c(4, 4))
for(i in 2:17){
  
  ylims=c(0,0.01)
  plot(density(orig_df[,i]), 
       main=paste0(colnames(orig_df)[i]),ylim=ylims)
}
dev.off()

thresh <- read.table("Distributions/thresholds/MLN_major_thresholds.txt", T, "\t")

orig_df$CD3[orig_df$CD3<thresh$CD3[1]] <- 0
orig_df$CD19[orig_df$CD19<thresh$CD19[1]] <- 0
orig_df$CD4[orig_df$CD4<thresh$CD4[1]] <- 0
orig_df$CD8[orig_df$CD8<thresh$CD8[1]] <- 0

#### write new df of CD4+ CD8+ CD19+ and CD3+
og_df <- read.table("inputs/MLN_df.txt", F, '\t')

#CD3
write.table(og_df[which(orig_df$CD3>thresh$CD3[1]),], "inputs/CD3_MLN_df.txt", 
            sep='\t', row.names=F, quote=F, col.names = F)

#CD19
write.table(og_df[which(orig_df$CD19>thresh$CD19[1]),], "inputs/CD19_MLN_df.txt", 
            sep='\t', row.names=F, quote=F, col.names = F)

#CD4
write.table(og_df[which(orig_df$CD4>thresh$CD4[1]),], "inputs/CD4_MLN_df.txt", 
            sep='\t', row.names=F, quote=F, col.names = F)

#CD8
write.table(og_df[which(orig_df$CD8>thresh$CD8[1]),], "inputs/CD8_MLN_df.txt", 
            sep='\t', row.names=F, quote=F, col.names = F)
