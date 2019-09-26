#### MLN Stimulation

### Clean up columns

### Evaluate for Plate Bias

setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/MLN_stimulations")

library(ggplot2)
library(reshape)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(gridExtra)

full_cyt <- read.table("CBA_MLN_Final_CD.txt", T, '\t')

metadata <- full_cyt[,1:6]

cyt_df <- full_cyt[,7:ncol(full_cyt)]
colnames(cyt_df) <- sub("A.\\.", "", colnames(cyt_df))
colnames(cyt_df) <- sub("A..\\.", "", colnames(cyt_df))
colnames(cyt_df) <- sub("B.\\.", "", colnames(cyt_df))
colnames(cyt_df) <- gsub("\\.\\.pg.ml.", "", colnames(cyt_df))
colnames(cyt_df) <- gsub("\\.pg.ml.", "", colnames(cyt_df))
colnames(cyt_df) <- gsub("\\.\\.", "_", colnames(cyt_df))
colnames(cyt_df) <- gsub("\\.", "_", colnames(cyt_df))

cyt_df <- apply(cyt_df, 2, function(x) gsub("<", "", x))
cyt_df <- apply(cyt_df, 2, function(x) gsub(">", "", x))
cyt_df <- apply(cyt_df, 2, function(x) gsub("ND", 0, x))
cyt_df <- apply(cyt_df, 2, function(x) gsub("--", 0, x))
cyt_df <- apply(cyt_df, 2, function(x) as.numeric(as.character(x)))

cyt_full <- cbind(metadata, cyt_df)

cyt_melt <- melt(cyt_full, measure.vars = colnames(cyt_full)[7:19], id.vars = colnames(cyt_full)[c(1:4,6)])
cyt_melt$value <- as.numeric(as.character(cyt_melt$value))
cyt_melt$challenge <- gsub("CD3/CD.*", "CD3_CD28", cyt_melt$challenge)
cyt_cast <- cast(cyt_melt, mouse_id ~ variable + challenge, mean)

cyt_cast$mouse_id <- gsub(" ", "", cyt_cast$mouse_id)
cyt_nums <- cyt_cast[,-1]
rownames(cyt_nums) = cyt_cast$mouse_id

write.table(cyt_cast, "MLN_stimulation_flat.txt", sep='\t', row.names=F, quote=F)
## now fix meta

meta <- cyt_full[,1:6]
meta$value <- sample(1:10000, nrow(meta))
meta <- aggregate(. ~ mouse_id + Location + Genetics + Plate, data=meta, mean)
meta$mouse_id <- gsub(" ", "", meta$mouse_id)
rownames(meta) <- meta$mouse_id
meta = meta[as.character(cyt_cast$mouse_id),1:4]

##

hh <- rowAnnotation(df = data.frame(Environment=meta$Location, Genetics=meta$Genetics, Plate=meta$Plate),
                    col = list(Environment = c(NYU = "mediumorchid3",
                                               SF = "red3"),
                               Genetics = c(AtgW = "darkorange1",
                                            AtgE = "dodgerblue2",
                                            AtgH = "navyblue",
                                            B6 = "darkorange1",
                                            NOD2 = "mediumspringgreen")))

df <- data.frame(names = colnames(cyt_nums))
df$og_names <- df$names
df$names <- gsub("_Bacillus", ":Bacillus", df$names)
df$names <- gsub("_Bacteroides", ":Bacteroides", df$names)
df$names <- gsub("_Candida", ":Candida", df$names)
df$names <- gsub("_CD3", ":CD3", df$names)
df$names <- gsub("_Clostri", ":Clostri", df$names)
df$names <- gsub("_PBS", ":PBS", df$names)
df$names <- gsub("_Pseudo", ":Pseudo", df$names)
df$names <- gsub("_Staph", ":Staph", df$names)
df$Cytokine <- gsub(":.*", "", df$names)
df$Challenge <- gsub(".*:", "", df$names)
ha <- columnAnnotation(df = data.frame(Cytokine=df$Cytokine))
labs = df$Challenge
ll <- columnAnnotation(labels = anno_text(labs, gp=gpar(fontsize=5),
                                          which = "column", rot=60, 
                                          just = 'right', offset=unit(1.2,"cm")), 
                       height = unit(1,"cm"))


pdf("cytokine_heat.pdf", height = 7, width = 11)
Heatmap(#t(scale(t(new_nums))),
  #scale(new_nums),
  scale(log2(cyt_nums+1)),
  #new_nums,
  show_row_names = T, show_column_names = F,
  top_annotation = ha,
  bottom_annotation = ll,
  bottom_annotation_height = unit(1.2,"cm"),
  heatmap_legend_param = list(title = "Cytokine\nLevel"),
  cluster_rows = T, cluster_columns = T, #row_names_side = 'left',
  row_names_gp = gpar(fontsize=6),
  #row_names_max_width = unit(10,'cm'),
  col = colorRamp2(c(0,10), c("white", "navyblue"))) + hh
  #col = colorRamp2(c(-2,0,2), c("blue", "white", "red"))) + hh
dev.off()

##

pc <- prcomp(log2(cyt_nums+1), scale = T)
summary(pc)
expl_var <- pc$sdev^2/sum(pc$sdev^2)*100

plot.data = data.frame(meta, cyt_nums, pc$x)
plot.data$Environment <- gsub("NYU", "lab", plot.data$Location)
plot.data$Environment <- gsub("SF", "wild", plot.data$Environment)
plot.data$Environment <- factor(plot.data$Environment, levels=c("lab","wild"))
plot.data$Genotype <- plot.data$Genetics
plot.data$Genotype <- gsub("AtgW", "B6", plot.data$Genotype)

g0 <- ggplot(plot.data, aes(PC1, PC2, color = Plate, shape = Genotype)) +
  geom_point(size = 3) + 
  ylab(paste0("PC2 ", round(expl_var[2],2), "% expl. variation")) +
  xlab(paste0("PC1 ", round(expl_var[1],2), "% expl. variation")) +
  #scale_color_manual(values = c("mediumpurple1", "red3")) +
  theme_bw()

pdf("cytokines_pca_plate.pdf", height = 10, width = 8)
grid.arrange(g0, g0 + facet_wrap(~ Genotype, ncol = 2), ncol = 1)
dev.off()

PC = pc
rownames(PC$x) <- 1:nrow(PC$x)
rownames(PC$rotation) <- colnames(cyt_nums)
data <- data.frame(obsnames=row.names(PC$x), PC$x[,1:2])
datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation[,1:2])
datapc$var1 <- scales::rescale(datapc$PC1, c(min(data$PC1),max(data$PC1)))
datapc$var2 <- scales::rescale(datapc$PC2, c(min(data$PC2),max(data$PC2)))

datapc$mult <- abs(datapc$PC1*datapc$PC2)
datapc <- datapc[order(datapc$mult, decreasing = T),]
datapc2 = datapc[1:20,]

g_seg1=ggplot(plot.data, aes(PC1,PC2, color = Environment, shape=Genotype, label=mouse_id)) +
  geom_point(size=3) + ggtitle("PCA of Cytokine Responses") +
  geom_text_repel(data=plot.data, aes(x=PC1, y=PC2, label=mouse_id), size = 3) +  
  ylab(paste0("PC2 ", round(expl_var[2],2), "% expl. variation")) +
  xlab(paste0("PC1 ", round(expl_var[1],2), "% expl. variation")) +
  scale_color_manual(values = c("mediumorchid3", "red3")) + coord_equal() +
  theme_bw()
g_seg=ggplot(plot.data, aes(PC1,PC2)) +
  geom_point(size=3, alpha=0) + ggtitle("Loadings for Cytokine+Stim") +
  scale_color_manual(values = c("mediumorchid3", "red3")) + coord_equal() + 
  geom_segment(data=datapc2, aes(x=0, y=0, xend=var1, yend=var2), 
               arrow=arrow(length=unit(0.2,"cm")), alpha=0.5) +
  geom_text_repel(data=datapc2, aes(x=var1, y=var2, label=varnames), size = 3) +
  theme_bw()
g_seg

pdf("cytokines_PCA_biplot.pdf", height = 12, width = 15)
grid.arrange(g_seg1+theme(legend.position = 'none'), g_seg, nrow = 2)
dev.off()
