##### ITS PCOA

setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/data/ITS")
### bray curtis

dist <- read.table("distance-matrix.tsv", T, '\t')
rownames(dist) <- dist$X
colnames(dist) <- dist$X
dist<-dist[,-1]

library(ape)
pcoa_obj <- pcoa(dist)
round(pcoa_obj$values$Relative_eig[1:2]*100,2)

plotter <- data.frame(PCoA_1=pcoa_obj$vectors[,1], PCoA_2=pcoa_obj$vectors[,2])

meta = read.table("metadata.tsv", T, '\t')
rownames(meta) <- meta$sample.id
meta=meta[rownames(plotter),]
plotter$Environment <- gsub(".*NYU.*", "lab", meta$SampleGroup)
plotter$Environment <- gsub(".*SF.*", "wild", plotter$Environment)
plotter$Environment <- factor(plotter$Environment, levels=c("lab","wild"))
plotter$Genotype <- gsub(".*AtgW.*", "WT", meta$Genotype)
plotter$Genotype <- gsub(".*B6.*", "B6", plotter$Genotype)
plotter$Genotype <- gsub("Atg16L1.T316A.T316A", "Atg16L1 T300A/T300A", plotter$Genotype)
#plotter$Genotype <- gsub("Atg16L1.T316A.WT.WT", "Atg16L1 T300A/+", plotter$Genotype)
plotter$Genotype <- gsub("Atg16L1.T316A.WT.WT", "B6", plotter$Genotype) ## apparently these are B6 animals?
plotter$Genotype <- gsub(".*NOD2.*", "NOD2 -/-", plotter$Genotype)

library(ggplot2)
g0 <- ggplot(plotter, aes(PCoA_1, PCoA_2, color = Environment))+ #, shape = Genotype)) +
  geom_point(size = 3) + 
  ylab(paste0("PCoA_2 ", "9.12", "% expl. variation")) +
  xlab(paste0("PCoA_1 ", "27.64", "% expl. variation")) +
  scale_color_manual(values = c("mediumpurple1", "red3")) +
  scale_shape_manual(values=c(17,15,3))+
  theme_bw()

library(gridExtra)
pdf("ITS_bray_pcoa_simp.pdf", height = 10, width = 8)
grid.arrange(g0, g0 + facet_wrap(~ Genotype, ncol = 2), ncol = 1)
dev.off()
