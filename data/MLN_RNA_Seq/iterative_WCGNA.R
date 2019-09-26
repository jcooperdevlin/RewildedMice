#### RNA_seq play and gender correction

library(matrixStats)
library(WGCNA)
library(dendextend)
library(fpc)
library(ggrepel)

setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/data/MLN_RNA_Seq")

counts <- read.table("normalizedCounts.txt", T, '\t')
samps <- gsub("X", "", colnames(counts)[-1])
gnames <- as.character(counts$GeneSymbol)
counts <- data.matrix(t(counts[,-1]))
colnames(counts)=gnames
rownames(counts)=samps
#counts <- counts[,order(colVars(counts), decreasing = T)][,1:10000]
counts <- counts[,colVars(counts)>1]

write.table(t(counts), "iWGCNA/input.txt", row.names=T, quote=F, sep='\t')

meta <- read.table("mice_metadata.11.19_mouse_id.txt", T, '\t')
rownames(meta) <- gsub("-", "_", as.character(meta$mouse_id))
meta<- meta[samps,]
meta$mouse_id <- rownames(meta)

nSets=1
multiExpr = vector(mode = "list", length = 1)
multiExpr[[1]] = list(data = as.data.frame(counts));
#names(multiExpr[[1]]$data) = femData$substanceBXH;
#rownames(multiExpr[[1]]$data) = names(femData)[-c(1:8)];

# Check that the data has the correct format for many functions operating on multiple sets:
exprSize = checkSets(multiExpr)

gsg = goodSamplesGenesMS(multiExpr, verbose = 3);


sampleTrees = list()
for (set in 1:nSets)
{
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}
plot(sampleTrees[[set]], main = paste("Sample clustering on all genes"),
     xlab="", sub="", cex = 0.7)

#cut by height of dendrogram for outliers
#cutHeights=list()
#cutHeights[[1]]=7000

#plot(sampleTrees[[set]], main = paste("Sample clustering on all genes"),
#     xlab="", sub="", cex = 0.7);
#abline(h=cutHeights[set], col = "red")


#for (set in 1:nSets)
#{
  # Find clusters cut by the line
#  labels = cutreeStatic(sampleTrees[[set]], cutHeight = cutHeights[set])
  # Keep the largest one (labeled by the number 1)
#  keep = (labels==1)
#  multiExpr[[set]]$data = multiExpr[[set]]$data[keep, ]
#}
#collectGarbage();
# Check the size of the leftover data
#exprSize = checkSets(multiExpr)
#exprSize


#

#
# meta

traitData = meta
dim(traitData)
names(traitData)
# remove columns that hold information we do not need.
allTraits = traitData[,c("mouse_id", "Environment", "Genotype")];
# See how big the traits are and what are the trait and sample names
dim(allTraits)
names(allTraits)

# Form a multi-set structure that will hold the clinical traits.
Traits = vector(mode="list", length = nSets);
for (set in 1:nSets)
{
  setSamples = rownames(multiExpr[[set]]$data);
  traitRows = match(setSamples, allTraits$mouse_id);
  Traits[[set]] = list(data = allTraits[traitRows, -1]);
  rownames(Traits[[set]]$data) = allTraits[traitRows, 1];
}
collectGarbage();
# Define data set dimensions
nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples;

#
##
setLabels = c("all")
shortLabels = c("all")
save(multiExpr, Traits, nGenes, nSamples, setLabels, shortLabels, exprSize,
     file = "Consensus-dataInput.RData")



#
#
# lets do this

lnames = load(file = "Consensus-dataInput.RData")
lnames
# Get the number of sets in the multiExpr structure.
nSets = checkSets(multiExpr)$nSets


# Choose a set of soft-thresholding powers
powers = c(seq(4,10,by=1), seq(12,20, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr2[[set]]$data, powerVector=powers,
                                                     verbose = 2)[[2]]);
collectGarbage();
# Plot the results:
colors = c("black", "red")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)
pdf(file = "scaleFreeAnalysis.pdf", wi = 8, he = 6)
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}
dev.off()

# we chose 6 based on this

net = blockwiseConsensusModules(
  multiExpr, power = 6, minModuleSize = 10, deepSplit = 4,
  maxBlockSize = 10000,
  pamRespectsDendro = FALSE,
  mergeCutHeight = 0.01, numericLabels = TRUE,
  minKMEtoStay = 0,
  saveTOMs = TRUE, verbose = 5)

consMEs = net$multiMEs;
moduleLabels = net$colors;
# Convert the numeric labels to color labels
moduleColors = labels2colors(moduleLabels)
consTree = net$dendrograms[[1]];
#A quick way to take a look at the results is to plot the gene dendrogram and the corresponding module colors:
sizeGrWindow(8,6);
pdf(file = "ConsensusDendrogram-auto.pdf", wi = 8, he = 6)
plotDendroAndColors(consTree, moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")
dev.off()




#
#

#### iterate through WGCNA to reduce gene sets into smaller sets

net = blockwiseConsensusModules(
  multiExpr, power = 6, minModuleSize = 10, deepSplit = 4,
  maxBlockSize = 11000,
  pamRespectsDendro = FALSE,
  mergeCutHeight = 0.001, numericLabels = TRUE,
  minKMEtoStay = 0,
  saveTOMs = TRUE, verbose = 4)

current_modules <- data.frame(Genes=colnames(multiExpr[[1]]$data), module = net$colors)
current_modules <- subset(current_modules, module != 0)

pammer <- function(gene_ids){
  newdata <- multiExpr[[1]]$data[,as.character(gene_ids$Genes)]
  
  new_dist=dist(t(newdata))
  new_clust = hclust(new_dist)
  
  new_pam <- pamk(new_dist)
  new_modules <- data.frame(Genes=colnames(newdata), module = new_pam$pamobject$clustering)
  new_modules$module <- paste0(mods[i], ".", new_modules$module)
  return(new_modules)
}

mods = unique(current_modules$module)
mod_max = max(mods)
updated_modules <- current_modules
for(i in 1:length(mods)){
  curr <- subset(updated_modules, module == mods[i])
  old <- subset(updated_modules, module != mods[i])
  if(nrow(curr)>100){
    it1=pammer(curr)
    updated_modules <- rbind(old, it1)
  }
}

mods = unique(updated_modules$module)
updated_modules2 <- current_modules
for(i in 1:length(mods)){
  curr <- subset(updated_modules2, module == mods[i])
  old <- subset(updated_modules2, module != mods[i])
  if(nrow(curr)>100){
    it1=pammer(curr)
    updated_modules2 <- rbind(old, it1)
  }
}

mods = unique(updated_modules2$module)
updated_modules3 <- updated_modules2
for(i in 1:length(mods)){
  curr <- subset(updated_modules3, module == mods[i])
  old <- subset(updated_modules3, module != mods[i])
  if(nrow(curr)>100){
    it1=pammer(curr)
    updated_modules3 <- rbind(old, it1)
  }
}

mods = unique(updated_modules3$module)
updated_modules4 <- updated_modules3
for(i in 1:length(mods)){
  curr <- subset(updated_modules4, module == mods[i])
  old <- subset(updated_modules4, module != mods[i])
  if(nrow(curr)>100){
    it1=pammer(curr)
    updated_modules4 <- rbind(old, it1)
  }
}

mods = unique(updated_modules4$module)
updated_modules5 <- updated_modules4
for(i in 1:length(mods)){
  curr <- subset(updated_modules5, module == mods[i])
  old <- subset(updated_modules5, module != mods[i])
  if(nrow(curr)>100){
    it1=pammer(curr)
    updated_modules5 <- rbind(old, it1)
  }
}

mods = unique(updated_modules5$module)
updated_modules6 <- updated_modules5
for(i in 1:length(mods)){
  curr <- subset(updated_modules6, module == mods[i])
  old <- subset(updated_modules6, module != mods[i])
  if(nrow(curr)>100){
    it1=pammer(curr)
    updated_modules6 <- rbind(old, it1)
  }
}

mods = unique(updated_modules6$module)
updated_modules7 <- updated_modules6
for(i in 1:length(mods)){
  curr <- subset(updated_modules7, module == mods[i])
  old <- subset(updated_modules7, module != mods[i])
  if(nrow(curr)>100){
    it1=pammer(curr)
    updated_modules7 <- rbind(old, it1)
  }
}

mods = unique(updated_modules7$module)
updated_modules8 <- updated_modules7
for(i in 1:length(mods)){
  curr <- subset(updated_modules8, module == mods[i])
  old <- subset(updated_modules8, module != mods[i])
  if(nrow(curr)>100){
    it1=pammer(curr)
    updated_modules8 <- rbind(old, it1)
  }
}

mods = unique(updated_modules8$module)
updated_modules9 <- updated_modules8
for(i in 1:length(mods)){
  curr <- subset(updated_modules9, module == mods[i])
  old <- subset(updated_modules9, module != mods[i])
  if(nrow(curr)>100){
    it1=pammer(curr)
    updated_modules9 <- rbind(old, it1)
  }
}

##
#
#
##### two heatmaps
##### one before module collapse with the module labels
##### and the second after


big_heat <- multiExpr[[1]]$data[,as.character(updated_modules9$Genes)]
small_heat <- aggregate(t(big_heat), by = list(updated_modules9$module), "mean")
colnames(small_heat)[1] <- "Module"
write.table(small_heat, "gene_modules_210.txt", sep="\t", row.names=F, quote=F)
write.table(updated_modules9, "gene_module_key.txt", sep='\t', row.names=F, quote=F)

big_heat = t(log2(big_heat+1))

png("Full_heat_9488_genes.png", height = 8, width = 10, units="in", res=300)
Heatmap(big_heat,
        col = colorRamp2(c(0,10), c("white", "purple4")),
        show_column_names = T, show_row_names = F)
dev.off()

pc <- prcomp(t(big_heat), scale=F)
plotter <- data.frame(PC1=pc$x[,1], PC2=pc$x[,2], id = colnames(big_heat))
plotter$env <- rep("red3", nrow(plotter))
plotter$env[grepl("_", plotter$id)] <- "mediumorchid3"

pc_before=ggplot(plotter, aes(PC1, PC2, label=id, color=env)) +
  geom_text_repel(show.legend=F) + 
  scale_color_manual(values=c("mediumorchid3", "red3")) +
  theme_bw()



png("small_heat_210_modules.png", height = 8, width = 10, units="in", res=300)
Heatmap(small_heat[,-1],
        col = colorRamp2(c(0,10), c("white", "purple4")),
        show_column_names = T, show_row_names = F)
dev.off()


pc <- prcomp(t(small_heat), scale=F)
plotter <- data.frame(PC1=pc$x[,1], PC2=pc$x[,2], id = colnames(small_heat))
plotter$env <- rep("red3", nrow(plotter))
plotter$env[grepl("_", plotter$id)] <- "mediumorchid3"

pc_after=ggplot(plotter, aes(PC1, PC2, label=id, color=env)) +
  geom_text_repel(show.legend=F) + 
  scale_color_manual(values=c("mediumorchid3", "red3")) +
  theme_bw()

pdf("pca_diff_modules.pdf", height=6, width = 12)
grid.arrange(pc_before, pc_after, nrow=1)
dev.off()

#
#
#
#
##
#

mods = unique(updated_modules$module)
mod_max = max(mods)
updated_modules2 <- updated_modules
for(i in 1:length(mods)){
  curr <- subset(updated_modules2, module == mods[i])
  old <- subset(updated_modules2, module != mods[i])
  if(nrow(curr)>100){
    newdata <- multiExpr[[1]]$data[,as.character(curr$Genes)]
    
    new_dist=dist(t(newdata))
    new_clust = hclust(new_dist)
    
    new_pam <- pamk(new_dist)
    print(paste0("Splitting by a k of ", new_pam$nc))
    
    new_modules <- data.frame(Genes=colnames(newdata), module = new_pam$pamobject$clustering)
    new_modules$module <- paste0(mods[i], ".", new_modules$module)
    
    updated_modules2 <- rbind(old, new_modules)
  }
}

table(current_modules$module)
table(updated_modules$module)
table(updated_modules2$module)





#
#
#

## are they enriched?

library(goseq)
library(GenomicFeatures)
library("Mus.musculus")
seqlengths(Mus.musculus)

DEgenes <- current_modules

genes <- rep(0, nrow(DEgenes))
genes[DEgenes$module==9] <- 1
names(genes) <- DEgenes$Gene


pwf=nullp(genes,"mm10","geneSymbol")
head(pwf)

GO.wall=goseq(pwf,"mm10","geneSymbol")
GO.wall$padj <- p.adjust(GO.wall$over_represented_pvalue, method ="BH")
GO.sig <- subset(GO.wall, over_represented_pvalue < 0.01)



