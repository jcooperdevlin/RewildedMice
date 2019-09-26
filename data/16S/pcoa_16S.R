##### 16S Figure

####

setwd("/Users/devlij03/Google Drive/Png/Ken_transfer/3.15/figs/3af/16S")

library(ggplot2)

pcoa <- read.table("pcoa_raw.txt", T, '\t')

pcoa$Environment <- gsub("NYU", "lab", pcoa$Location)
pcoa$Environment <- gsub("SF", "wild", pcoa$Environment)
pcoa$Environment <- factor(pcoa$Environment, levels=c("lab","wild"))
pcoa$Genotype <- pcoa$Genetics
pcoa$Genotype <- gsub("AtgW", "B6", pcoa$Genotype)

g0 <- ggplot(pcoa, aes(PCoA_1, PCoA_2, color = Environment, shape = Genotype)) +
  geom_point(size = 3) + 
  ylab(paste0("PC2 ", "17.71", "% expl. variation")) +
  xlab(paste0("PC1 ", "24.61", "% expl. variation")) +
  scale_color_manual(values = c("mediumpurple1", "red3")) +
  theme_bw()


library(plotly)
p <- plot_ly(pcoa, x = ~PCoA_1, y = ~PCoA_2, z = ~PCoA_3, 
             color = pcoa$Environment,
             type = 'scatter3d', mode = 'markers', text = pcoa$Environment,
             hoverinfo = "text") %>% layout(paper_bgcolor='transparent', font = t)
p

library(scatterplot3d)
PC1 = pcoa$PCoA_1
PC2 = pcoa$PCoA_2
PC3 = pcoa$PCoA_3

pcoa$color <- gsub("wild", "red3", pcoa$Environment)
pcoa$color <- gsub("lab", "mediumorchid3", pcoa$color)

pcoa$symbol <- gsub("AtgE", 16, pcoa$Genotype)
pcoa$symbol <- gsub("AtgH", 17, pcoa$symbol)
pcoa$symbol <- gsub("B6", 15, pcoa$symbol)
pcoa$symbol <- gsub("NOD2", 3, pcoa$symbol)
pcoa$symbol <- as.numeric(pcoa$symbol)

pdf("weighted_unifrac_pcoa.pdf", height = 7, width = 8)
scatterplot3d(x =PC3, y=PC2, z=PC1,
              pch = pcoa$symbol, cex.symbols = 1.5, box=F, 
              color = pcoa$color)
dev.off()




#
#
#
### bray curtis

dist <- read.table("bray_curtis_distance-matrix.tsv", T, '\t')
rownames(dist) <- dist$X
colnames(dist) <- dist$X
dist<-dist[,-1]

library(ape)
pcoa_obj <- pcoa(dist)
round(pcoa_obj$values$Relative_eig[1:2]*100,2)

plotter <- data.frame(PCoA_1=pcoa_obj$vectors[,1], PCoA_2=pcoa_obj$vectors[,2])

plotter$Environment <- gsub(".*NYU.*", "lab", rownames(dist))
plotter$Environment <- gsub(".*SF.*", "wild", plotter$Environment)
plotter$Environment <- factor(plotter$Environment, levels=c("lab","wild"))
plotter$Genotype <- gsub(".*AtgW.*", "B6", rownames(dist))
plotter$Genotype <- gsub(".*B6.*", "B6", plotter$Genotype)
plotter$Genotype <- gsub(".*AtgE.*", "AtgE", plotter$Genotype)
plotter$Genotype <- gsub(".*AtgH.*", "AtgH", plotter$Genotype)
plotter$Genotype <- gsub(".*NOD2.*", "NOD2", plotter$Genotype)
                         

g0 <- ggplot(plotter, aes(PCoA_1, PCoA_2, color = Environment, shape = Genotype)) +
  geom_point(size = 3) + 
  ylab(paste0("PCoA_2 ", "12.12", "% expl. variation")) +
  xlab(paste0("PCoA_1 ", "19.08", "% expl. variation")) +
  scale_color_manual(values = c("mediumpurple1", "red3")) +
  theme_bw()

pdf("16S_bray_pcoa.pdf", height = 10, width = 8)
grid.arrange(g0, g0 + facet_wrap(~ Genotype, ncol = 2), ncol = 1)
dev.off()
