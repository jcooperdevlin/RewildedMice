#### RNA_seq play and gender correction

library(matrixStats)
library(WGCNA)

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
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
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
