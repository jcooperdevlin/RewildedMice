### serum radar 6.27
library(reshape)
library(ggplot2)

setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/data/Serum_stimulations")

serum_df <- read.table("Serum_names_final.txt", T)
meta <- read.table("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/data/metadata/mice_metadata.11.19_mouse_id.txt", 
                   header=T, sep='\t')
rownames(meta) <- as.character(meta$mouse_id)
cyt_keep <- intersect(unique(serum_df$sample), meta$mouse_id)
meta <- meta[cyt_keep,]
serum_df$Environment <- meta$Environment
serum_df$Genotype <- meta$Genotype

ser_melt <- melt(serum_df, measure.vars = colnames(serum_df)[2:11], id.vars = colnames(serum_df)[c(1,12,13)])
plotter_slim <- aggregate(ser_melt$value, list(ser_melt$Environment, ser_melt$Genotype, ser_melt$variable),mean)
colnames(plotter_slim) <- c("Environment", "Genotype", "variable", "value")

plotter_slim$Environment <- factor(plotter_slim$Environment, levels = c("wild", "lab"))

g=ggplot(plotter_slim, aes(x=variable, y=log2(value+1))) +
  geom_polygon(aes(group = Environment, color = Environment,fill = Environment),
               alpha=0.4, size = 1, show.legend = TRUE) +
  coord_polar(clip='on')+
  #scale_y_continuous(limits = c(0,6)) +
  #geom_text(data=plotter_slim, aes(x=cyto, y=labelery, label=labelerx), size = 5) + 
  scale_color_manual(values = c("red3", "mediumpurple1"))+
  scale_fill_manual(values = c("red3", "mediumpurple1"))+
  xlab("Cytokines") + ylab("Log2\nCytokine Level") +
  theme(
    panel.background = element_rect(fill='white'),
    #plot.margin=margin(10, 10, 5, 5),
    panel.grid.major = element_line(color='grey92'),
    #axis.text.x=element_blank(),
    axis.text.y = element_text(size = 12, color = 'black'),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    legend.title = element_text(size=15))

pdf("serum_diff.pdf", height = 7, width = 8)
g+facet_wrap(~Genotype)
dev.off()



######
source("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/DIST_EffectSize/var_plotter.R")

serum_df <- read.table("Serum_names_final.txt", T)
meta <- read.table("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/data/metadata/mice_metadata.11.19_mouse_id.txt", 
                   header=T, sep='\t')
rownames(meta) <- as.character(meta$mouse_id)
cyt_keep <- intersect(unique(serum_df$sample), meta$mouse_id)
meta <- meta[cyt_keep,]

input=log2(serum_df[,-1]+1)
rownames(input) <- serum_df[,1]

rowZeros <- apply(input, 1, function(x) {sum(x == 0)})
colZeros <- apply(input, 2, function(x) {sum(x == 0)})

input = input[which(rowZeros<0.5*ncol(input)),which(colZeros<0.5*nrow(input)) ]

meta_keep<-meta[rownames(input),]

effectors=meta_keep[,c(2,3,4,6,8)]
col_var=meta_keep$Environment
shape_var=meta_keep$Genotype

gg_ser_cyt <- var_plotter(input, effectors, col_var, "Environment",
                          shape_var, "Genotype", "Serum Cytokines")

pdf("serum_cytokine_pcoa_ES.pdf", height = 5, width = 15)
grid.arrange(gg_ser_cyt)
dev.off()
