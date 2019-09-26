#### Cytokine Stimulation Correlation and PCA

### 3.15.19 update meeting

library(ComplexHeatmap)
library(psych)
library(reshape)
library(circlize)
library(ggplot2)
library(gridExtra)
library(ggrepel)

string.to.colors = function(string, colors=NULL){
  if (is.factor(string)){
    string = as.character(string)
  }
  if (!is.null(colors)){
    #if (length(colors)!=length(unique(string))){
    #  break("The number of colors must be equal to the number of unique elements.")
    #}
    #else {
    conv = cbind(unique(string), colors)
    #}
  } else {
    conv = cbind(unique(string), rainbow(length(unique(string))))
  }
  unlist(lapply(string, FUN=function(x){conv[which(conv[,1]==x),2]}))
}


setwd("/Users/devlij03/Google Drive/Png/Ken_transfer/3.15/figs/2a")


stim_df <- read.table("MLN_stimulations_names.txt", header=T, sep='\t')
stim_df$challenge <- gsub("CD3/CD.*", "CD3/CD28", stim_df$challenge)

stim_melt <- melt(stim_df, measure.vars = colnames(stim_df)[7:19], id.vars = colnames(stim_df)[1:4])
stim_melt$value <- gsub("<.*", 0, stim_melt$value)
stim_melt$value <- as.numeric(stim_melt$value)
stim_melt$value[is.na(stim_melt$value)] <- 0
stim_cast <- cast(stim_melt, mouse_id ~ variable + challenge, mean)

stim_cast$mouse_id <- gsub(" ", "", stim_cast$mouse_id)

write.table(stim_cast, "MLN_stimulation_flat.txt", row.names = F, sep='\t', quote=F)
#nums <- apply(stim_df[,7:19], 2, function(x) gsub("<.*", 0, x))
#nums <- apply(nums, 2, as.numeric)
#nums[is.na(nums)] <- 0

nums <- stim_cast[,2:ncol(stim_cast)]
cc <- colnames(nums)
rr <- rownames(nums)
#nums <- scale(nums)
nums <- data.frame(data.matrix(nums))
colnames(nums) <- cc
rownames(nums) <- as.character(stim_cast$mouse_id)

meta <- stim_df[,1:6]
meta$value <- sample(1:10000, nrow(meta))
meta <- aggregate(. ~ mouse_id + Location + Genetics + Plate, data=meta, mean)
rownames(meta) <- gsub(" ", "", as.character(meta$mouse_id))
meta$mouse_id <- gsub(" ", "", as.character(meta$mouse_id))

meta = meta[as.character(stim_cast$mouse_id),1:4]

table(meta$Genetics, meta$Plate)

#new_nums <- matrix(nums[1,], nrow = 1)
#colnames(new_nums) <- colnames(nums)
#plates <- unique(stim_df$Plate)
#for (i in 1:length(plates)){
#  p1 <- nums[which(stim_df$Plate == plates[i]),]
#  p2 <- scale(p1)
#  new_nums <- rbind(new_nums, p2)
#}
#new_nums <- new_nums[-1,]
#new_nums[is.nan(new_nums)] <-0

new_nums <- nums

new_nums2 <- data.frame(matrix(rep(0,104), 1,104))
colnames(new_nums2) <- colnames(new_nums)
plates <- levels(meta$Plate)
for (i in 1:length(plates)){
  use <- new_nums[which(meta$Plate == plates[i]),]
  #use <- t(scale(t(use)))
  #use <- scale(use)
  #use <- log2(use+1)
  new_nums2<- rbind(new_nums2, use)
}
new_nums2<-new_nums2[-1,]
new_nums=new_nums2[as.character(meta$mouse_id),]
new_nums[is.na(new_nums)]<-0

##pause for heat

hh <- rowAnnotation(df = data.frame(Environment=meta$Location, Genetics=meta$Genetics, Plate=meta$Plate),
                    col = list(Environment = c(NYU = "mediumorchid3",
                                               SF = "red3"),
                               Genetics = c(AtgW = "darkorange1",
                                            AtgE = "dodgerblue2",
                                            AtgH = "navyblue",
                                            B6 = "darkorange1",
                                            NOD2 = "mediumspringgreen")))

df <- data.frame(names = colnames(new_nums))
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


pdf("cytokine_heat2.pdf", height = 5, width = 11)
Heatmap(#t(scale(t(new_nums))),
        #scale(new_nums),
        log2(new_nums+1),
        #new_nums,
        show_row_names = F, show_column_names = F,
        top_annotation = ha,
        bottom_annotation = ll,
        bottom_annotation_height = unit(1.2,"cm"),
        heatmap_legend_param = list(title = "Cytokine\nLevel"),
        cluster_rows = T, cluster_columns = T, #row_names_side = 'left',
        #row_names_gp = gpar(fontsize=8),
        #row_names_max_width = unit(10,'cm'),
        col = colorRamp2(c(0,10), c("white", "navyblue"))) + hh
        #col = colorRamp2(c(-2,0,2), c("blue", "white", "red"))) + hh
dev.off()

##


rownames(new_nums) <- 1:nrow(new_nums)
pc <- prcomp(t(scale(t(new_nums))), scale = F)
#pc <- prcomp(log2(new_nums+1), scale = F)
#pc <- prcomp(new_nums, scale = T)
summary(pc)
expl_var <- pc$sdev^2/sum(pc$sdev^2)*100


plot.data = data.frame(meta, new_nums, pc$x)
plot.data$Environment <- gsub("NYU", "lab", plot.data$Location)
plot.data$Environment <- gsub("SF", "wild", plot.data$Environment)
plot.data$Environment <- factor(plot.data$Environment, levels=c("lab","wild"))
plot.data$Genotype <- plot.data$Genetics
plot.data$Genotype <- gsub("AtgW", "B6", plot.data$Genotype)


g0 <- ggplot(plot.data, aes(PC1, PC2, color = Environment, shape = Genotype)) +
  geom_point(size = 6) + 
  ylab(paste0("PC2 ", round(expl_var[2],2), "% expl. variation")) +
  xlab(paste0("PC1 ", round(expl_var[1],2), "% expl. variation")) +
  scale_color_manual(values = c("mediumpurple1", "red3")) +
  theme_bw()+
  theme(axis.text = element_text(size=12, color="black"),
        axis.title = element_text(size=13, color="black"),
        legend.text = element_text(size=12, color="black"),
        legend.title = element_text(size=13, color="black"),
        plot.title = element_text(size=14, face='bold'))

pdf("cytokines_pca.pdf", height = 10, width = 8)
grid.arrange(g0, g0 + facet_wrap(~ Genotype, ncol = 2), ncol = 1)
dev.off()


PC = pc
rownames(PC$x) <- 1:nrow(PC$x)
rownames(PC$rotation) <- colnames(new_nums)
data <- data.frame(obsnames=row.names(PC$x), PC$x[,1:2])
datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation[,1:2])
datapc$var1 <- scales::rescale(datapc$PC1, c(min(data$PC1),max(data$PC1)))
datapc$var2 <- scales::rescale(datapc$PC2, c(min(data$PC2),max(data$PC2)))

datapc$mult <- abs(datapc$PC1*datapc$PC2)
datapc <- datapc[order(datapc$mult, decreasing = T),]
datapc2 = datapc[1:25,]

g_seg1=ggplot(plot.data, aes(PC1,PC2, color = Environment, shape=Genotype)) +
  geom_point(size=3) + ggtitle("PCA of Cytokine Responses") +
  ylab(paste0("PC2 ", round(expl_var[2],2), "% expl. variation")) +
  xlab(paste0("PC1 ", round(expl_var[1],2), "% expl. variation")) +
  scale_color_manual(values = c("mediumorchid3", "red3")) + coord_equal() +
  theme_bw()+
  theme(axis.text = element_text(size=12, color="black"),
        axis.title = element_text(size=13, color="black"),
        legend.text = element_text(size=12, color="black"),
        legend.title = element_text(size=13, color="black"),
        plot.title = element_text(size=14, face='bold'))
g_seg=ggplot(plot.data, aes(PC1,PC2)) +
  geom_point(size=3, alpha=0) + ggtitle("Loadings for Cytokine+Stim") +
  scale_color_manual(values = c("mediumorchid3", "red3")) + coord_equal() + 
  geom_segment(data=datapc2, aes(x=0, y=0, xend=var1, yend=var2), 
               arrow=arrow(length=unit(0.2,"cm")), alpha=0.5) +
  geom_text_repel(data=datapc2, aes(x=var1, y=var2, label=varnames), size = 4) +
  theme_bw()+
  theme(axis.text = element_text(size=12, color="black"),
        axis.title = element_text(size=13, color="black"),
        legend.text = element_text(size=12, color="black"),
        legend.title = element_text(size=13, color="black"),
        plot.title = element_text(size=14, face='bold'))
g_seg


pdf("cytokines_PCA_biplot.pdf", height = 12, width = 10)
grid.arrange(g_seg1+theme(legend.position = 'none'), g_seg, nrow = 2)
dev.off()


#
### pause for effect size along pc1
effector <- function(df){
  es <- c()
  type <- c()
  for (i in 2:c(ncol(df)-1)){
    if (is.character(df[,i]) | is.factor(df[,i])){
      if(length(unique(df[,i]))>1){
        var_form <- as.formula(paste("index~", colnames(df)[i]))
        fit <- unlist(round(sqrt(eta_sq(aov(var_form, df))[2]), 3))
        es <- c(es, fit)
        type <- c(type, colnames(df)[i])
      } else {
        es <- c(es, NA)
        type <- c(type, colnames(df)[i])
      }
    }
    if (is.numeric(df[,i])){
      fit <- round(cor(df$index, df[,i], method = "spearman"), 3)
      es <- c(es, fit)
      if(colnames(df)[i]=="Weight_gain."){
        type <- c(type, "WeightGain")
      } else {
        type <- c(type, "CytStim")
      }
    }
  }
  res <- data.frame(Variable = colnames(df[,2:c(ncol(df)-1)]),
                    EffectSize = es, type = type)
  return(res)
}

type_cols <- c(Plate="red2", Location="dodgerblue2", CytStim = "forestgreen",
               Wedge_cage="navyblue",Flow.date="grey50",
               Gender="hotpink",Genetics="darkorange1", Pregnant="purple", WeightGain="black", Genotype="darkorange1")

##
effect_check <- plot.data[,c(1:108,110)]
#effect_check <- plot.data[,c(1:109)]
colnames(effect_check)[109] <- "index"
#colnames(effect_check)[6] <- "Wedge_cage"
#effect_check$Wedge_cage <- factor(effect_check$Wedge_cage)

es_df <- effector(effect_check)

es_df <- es_df[order(es_df$EffectSize, decreasing = T),]
es_df$Variable <- factor(es_df$Variable, levels = rev(es_df$Variable))
es_df <- subset(es_df, !is.na(EffectSize))
es_df <- es_df[order(abs(es_df$EffectSize), decreasing = T),]
es_df$Variable <- factor(es_df$Variable, levels = rev(es_df$Variable))

es_filt <- subset(es_df, type != "Plate" & type != "CytStim")
es_cyt <- subset(es_df, type == "CytStim" & abs(EffectSize)>0.2)

es_filt$type <- gsub("Genetics", "Genotype", es_filt$type)

es_df <- rbind(es_filt, es_cyt)
gg_combo=ggplot(es_df, aes(Variable, abs(EffectSize), fill = type)) +
  geom_col() +
  guides(fill=guide_legend(title="Feature"))+
  ggtitle("EffectSize Along PC2")+
  ylab("EffectSize")+
  scale_fill_manual(values=type_cols)+
  theme_bw()+
  coord_flip() +
  theme(axis.text = element_text(size=12, color="black"),
        axis.title = element_text(size=13, color="black"),
        legend.text = element_text(size=12, color="black"),
        legend.title = element_text(size=13, color="black"),
        plot.title = element_text(size=14, face='bold'))

pdf("effect_size_calc_all_filt.pdf", height = 4, width = 6)
gg_combo
dev.off()


##lab
effect_check <- subset(plot.data[,c(1:108,110)], Location=="NYU")
#effect_check <- subset(plot.data[,c(1:109)], Location=="NYU")
colnames(effect_check)[109] <- "index"
#colnames(effect_check)[6] <- "Wedge_cage"
#effect_check$Wedge_cage <- factor(effect_check$Wedge_cage)


es_df <- effector(effect_check)

es_df <- es_df[order(es_df$EffectSize, decreasing = T),]
es_df$Variable <- factor(es_df$Variable, levels = rev(es_df$Variable))

es_df <- subset(es_df, !is.na(EffectSize))
gg_lab=ggplot(es_df, aes(Variable, EffectSize, fill = type)) +
  geom_col() +
  ggtitle("Lab Mice")+
  scale_fill_manual(values=type_cols)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=unit(1,"cm")))

###
##wild
effect_check <- subset(plot.data[,c(1:108,110)], Location=="SF")
#effect_check <- subset(plot.data[,c(1:109)], Location=="SF")
colnames(effect_check)[109] <- "index"
#colnames(effect_check)[6] <- "Wedge_cage"
#effect_check$Wedge_cage <- factor(effect_check$Wedge_cage)


es_df <- effector(effect_check)

es_df <- es_df[order(es_df$EffectSize, decreasing = T),]
es_df$Variable <- factor(es_df$Variable, levels = rev(es_df$Variable))

es_df <- subset(es_df, !is.na(EffectSize))

gg_wild=ggplot(es_df, aes(Variable, EffectSize, fill = type)) +
  geom_col() +
  ggtitle("Wild Mice")+
  scale_fill_manual(values=type_cols)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=unit(1,"cm")))


pdf("effect_size_calc_cyt2.pdf", height = 12, width = 14)
grid.arrange(gg_combo, gg_lab, gg_wild, nrow=3)
dev.off()

#
#
###lab pca
lab_nums<-new_nums[which(meta$Location=="NYU"),]
lab_meta <- meta[which(meta$Location=="NYU"),]
rownames(lab_nums) <- 1:nrow(lab_nums)
pc <- prcomp(t(scale(t(lab_nums))), scale = F)
summary(pc)
expl_var <- pc$sdev^2/sum(pc$sdev^2)*100


plot.data = data.frame(lab_meta, lab_nums, pc$x)
plot.data$Environment <- gsub("NYU", "lab", plot.data$Location)
plot.data$Environment <- gsub("SF", "wild", plot.data$Environment)
plot.data$Environment <- factor(plot.data$Environment, levels=c("lab","wild"))
plot.data$Genotype <- plot.data$Genetics
plot.data$Genotype <- gsub("AtgW", "B6", plot.data$Genotype)


g0 <- ggplot(plot.data, aes(PC1, PC2, color = Environment, shape = Genotype)) +
  geom_point(size = 3) + 
  ylab(paste0("PC2 ", round(expl_var[2],2), "% expl. variation")) +
  xlab(paste0("PC1 ", round(expl_var[1],2), "% expl. variation")) +
  scale_color_manual(values = c("mediumpurple1", "red3")) +
  theme_bw()

pdf("lab_cytokines_pca.pdf", height = 10, width = 8)
grid.arrange(g0, g0 + facet_wrap(~ Genotype, ncol = 2), ncol = 1)
dev.off()


PC = pc
rownames(PC$x) <- 1:nrow(PC$x)
rownames(PC$rotation) <- colnames(lab_nums)
data <- data.frame(obsnames=row.names(PC$x), PC$x[,1:2])
datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation[,1:2])
datapc$var1 <- scales::rescale(datapc$PC1, c(min(data$PC1),max(data$PC1)))
datapc$var2 <- scales::rescale(datapc$PC2, c(min(data$PC2),max(data$PC2)))

datapc$mult <- abs(datapc$PC1*datapc$PC2)
datapc <- datapc[order(datapc$mult, decreasing = T),]
datapc2 = datapc[1:20,]

lab_seg1=ggplot(plot.data, aes(PC1,PC2, color = Environment, shape=Genotype)) +
  geom_point(size=3) + ggtitle("PCA of Cytokine Responses") +
  ylab(paste0("PC2 ", round(expl_var[2],2), "% expl. variation")) +
  xlab(paste0("PC1 ", round(expl_var[1],2), "% expl. variation")) +
  scale_color_manual(values = c("mediumorchid3", "red3")) + coord_equal() +
  #ylim(c(-10,20))+xlim(c(-25,10))+
  theme_bw()
lab_seg=ggplot(plot.data, aes(PC1,PC2)) +
  geom_point(size=3, alpha=0) + ggtitle("Loadings for Cytokine+Stim") +
  scale_color_manual(values = c("mediumorchid3", "red3")) + coord_equal() + 
  geom_segment(data=datapc2, aes(x=0, y=0, xend=var1, yend=var2), 
               arrow=arrow(length=unit(0.2,"cm")), alpha=0.5) +
  geom_text_repel(data=datapc2, aes(x=var1, y=var2, label=varnames), size = 3) +
  #ylim(c(-10,20))+xlim(c(-25,10))+
  theme_bw()
lab_seg


#pdf("lab_cytokines_PCA_biplot.pdf", height = 12, width = 15)
#grid.arrange(g_seg1+theme(legend.position = 'none'), g_seg, nrow = 2)
#dev.off()



#
#
#
#
###wild pca
wild_nums<-new_nums[which(meta$Location=="SF"),]
wild_meta <- meta[which(meta$Location=="SF"),]
rownames(wild_nums) <- 1:nrow(wild_nums)
pc <- prcomp(t(scale(t(wild_nums))), scale = F)
summary(pc)
expl_var <- pc$sdev^2/sum(pc$sdev^2)*100


plot.data = data.frame(wild_meta, wild_nums, pc$x)
plot.data$Environment <- gsub("NYU", "lab", plot.data$Location)
plot.data$Environment <- gsub("SF", "wild", plot.data$Environment)
plot.data$Environment <- factor(plot.data$Environment, levels=c("lab","wild"))
plot.data$Genotype <- plot.data$Genetics
plot.data$Genotype <- gsub("AtgW", "B6", plot.data$Genotype)


g0 <- ggplot(plot.data, aes(PC1, PC2, color = Environment, shape = Genotype)) +
  geom_point(size = 3) + 
  ylab(paste0("PC2 ", round(expl_var[2],2), "% expl. variation")) +
  xlab(paste0("PC1 ", round(expl_var[1],2), "% expl. variation")) +
  scale_color_manual(values = c("red3")) +
  theme_bw()

pdf("wild_cytokines_pca.pdf", height = 10, width = 8)
grid.arrange(g0, g0 + facet_wrap(~ Genotype, ncol = 2), ncol = 1)
dev.off()


PC = pc
rownames(PC$x) <- 1:nrow(PC$x)
rownames(PC$rotation) <- colnames(wild_nums)
data <- data.frame(obsnames=row.names(PC$x), PC$x[,1:2])
datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation[,1:2])
datapc$var1 <- scales::rescale(datapc$PC1, c(min(data$PC1),max(data$PC1)))
datapc$var2 <- scales::rescale(datapc$PC2, c(min(data$PC2),max(data$PC2)))

datapc$mult <- abs(datapc$PC1*datapc$PC2)
datapc <- datapc[order(datapc$mult, decreasing = T),]
datapc2 = datapc[1:20,]

wild_seg1=ggplot(plot.data, aes(PC1,PC2, color = Environment, shape=Genotype)) +
  geom_point(size=3) + ggtitle("PCA of Cytokine Responses") +
  ylab(paste0("PC2 ", round(expl_var[2],2), "% expl. variation")) +
  xlab(paste0("PC1 ", round(expl_var[1],2), "% expl. variation")) +
  scale_color_manual(values = c("red3")) + coord_equal() +
  #ylim(c(-10,20))+xlim(c(-25,10))+
  theme_bw()
wild_seg=ggplot(plot.data, aes(PC1,PC2)) +
  geom_point(size=3, alpha=0) + ggtitle("Loadings for Cytokine+Stim") +
  scale_color_manual(values = c("mediumorchid3", "red3")) + coord_equal() + 
  geom_segment(data=datapc2, aes(x=0, y=0, xend=var1, yend=var2), 
               arrow=arrow(length=unit(0.2,"cm")), alpha=0.5) +
  geom_text_repel(data=datapc2, aes(x=var1, y=var2, label=varnames), size = 3) +
  #ylim(c(-10,20))+xlim(c(-25,10))+
  theme_bw()
wild_seg


#pdf("wild_cytokines_PCA_biplot.pdf", height = 12, width = 15)
#grid.arrange(g_seg1+theme(legend.position = 'none'), g_seg, nrow = 2)
#dev.off()

pdf("lab_wild_cytokines_PCA_biplot.pdf", height = 12, width = 15)
grid.arrange(lab_seg1+theme(legend.position = 'none'),
             wild_seg1+theme(legend.position = 'none'), 
             lab_seg, wild_seg, nrow = 2)
dev.off()

#
#
#
#
#
#
#
#### Some correlations

corr_heatmapper <- function(input, filename="corr.pdf"){
hah2 = corr.test(input, method = "spearman")
corr_mat = as.matrix(hah2$r)

orig_p = hah2$p
orig_p[col(orig_p) == row(orig_p) | upper.tri(orig_p)] <- NA
new_ps = p.adjust(orig_p, method = "BH", n = length(orig_p))
new_p_mat = matrix(new_ps, ncol = ncol(orig_p), nrow = nrow(orig_p))

new_p_mat[upper.tri(new_p_mat)] = t(new_p_mat)[upper.tri(new_p_mat)]
diag(new_p_mat) = rep(1,nrow(new_p_mat))

p_list = c(new_p_mat)
p_dims = dim(new_p_mat)

sym_list = c()
for (i in p_list){
  if (is.na(i)){
    sym = ""
  }
  else {
    if (i > 0){
      sym = ""
      if (i < 0.05){
        sym = "*"
        if (i<0.01){
          sym = "**"
          if(i<0.001){
            sym = "***"
          }
        }
      }
    }
  }
  sym_list = c(sym_list, sym)
}

sym_mat = matrix(sym_list, ncol = p_dims[1], nrow = p_dims[2])
colnames(sym_mat) <- colnames(corr_mat)
rownames(sym_mat) <- rownames(corr_mat)

corr_mat[is.na(corr_mat)] <- 0

corr_mat2=corr_mat[hclust(dist(corr_mat))$labels[hclust(dist(corr_mat))$order],
                   hclust(dist(corr_mat))$labels[hclust(dist(corr_mat))$order]]

namer2 <- gsub("corr", "dend", filename)
pdf(namer2, height = 7, width = 10)
plot(as.dendrogram(hclust(dist(corr_mat))))
dev.off()

df <- data.frame(names = colnames(corr_mat2))
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

#df <- df[order(df$Challenge, df$Cytokine),]
#corr_mat2 <- corr_mat[as.character(df$og_names),as.character(df$og_names)]
corr_mat2[upper.tri(corr_mat2, diag = F)] = 0

va <- rowAnnotation(df = df[,3:4],
                    col = list(Cytokine=c(IL.1a = "red", IL.5 = "blue", IFN.y = "green",
                                          TNF.a = "pink", KC_CXCL1 = "purple", MCP.1_CCL2 = "orange",
                                          IL.1b = "yellow", IL.10 = "turquoise",
                                          MIP.1a_CCL3 = "grey", MIP.1b_CCL4 = "black", IL.17 = "beige",
                                          IL.6 = "brown", IL.13 = "green4"),
                               Challenge=c(Bacillus.S.="red", Bacteroides.V. = "blue",
                                           Candida.A.="yellow", CD3.CD28="green",
                                           Clostridium.P.="pink", PBS="purple", Pseudomonas.A.="orange",
                                           Staph.A.="black")))

ha <- columnAnnotation(df = df[,3:4],
                    col = list(Cytokine=c(IL.1a = "red", IL.5 = "blue", IFN.y = "green",
                                          TNF.a = "pink", KC_CXCL1 = "purple", MCP.1_CCL2 = "orange",
                                          IL.1b = "yellow", IL.10 = "turquoise",
                                          MIP.1a_CCL3 = "grey", MIP.1b_CCL4 = "black", IL.17 = "beige",
                                          IL.6 = "brown", IL.13 = "green4"),
                               Challenge=c(Bacillus.S.="red", Bacteroides.V. = "blue",
                                           Candida.A.="yellow", CD3.CD28="green",
                                           Clostridium.P.="pink", PBS="purple", Pseudomonas.A.="orange",
                                           Staph.A.="black")))

hah <- Heatmap(corr_mat2, 
               col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
               top_annotation = ha,
               #split = df$Challenge,
               cluster_columns = F, cluster_rows = F,
               show_row_names = T, show_column_names = F,
               row_names_gp = gpar(fontsize = 6),
               row_names_side = "left",
               heatmap_legend_param=list(title = "Correlation"))+va

namer = filename
pdf(namer, height = 8, width = 8)
draw(hah)
dev.off()
}

nums_2_use <- plot.data[,4:107]
corr_heatmapper(nums_2_use, filename = "cyt_corr_all.pdf")

# nums_2_use <- subset(plot.data[,4:107], plot.data$Environment == "wild")
# corr_heatmapper(nums_2_use, filename = "cyt_corr_wild.pdf")
# 
# nums_2_use <- subset(plot.data[,4:107], plot.data$Environment == "lab")
# corr_heatmapper(nums_2_use, filename = "cyt_corr_lab.pdf")
# 
# nums_2_use <- subset(plot.data[,4:107], plot.data$Environment == "wild" & plot.data$Genotype == "B6")
# corr_heatmapper(nums_2_use, filename = "cyt_corr_wildB6.pdf")
# 
# nums_2_use <- subset(plot.data[,4:107], plot.data$Environment == "lab" & plot.data$Genotype == "B6")
# corr_heatmapper(nums_2_use, filename = "cyt_corr_labB6.pdf")
# 
# 
# nums_2_use <- subset(plot.data[,4:107], plot.data$Environment == "wild" & plot.data$Genotype == "AtgH")
# corr_heatmapper(nums_2_use, filename = "cyt_corr_wildAtgH.pdf")
# 
# nums_2_use <- subset(plot.data[,4:107], plot.data$Environment == "lab" & plot.data$Genotype == "AtgH")
# corr_heatmapper(nums_2_use, filename = "cyt_corr_labAtgH.pdf")
# 
# 
# nums_2_use <- subset(plot.data[,4:107], plot.data$Environment == "wild" & plot.data$Genotype == "AtgE")
# corr_heatmapper(nums_2_use, filename = "cyt_corr_wildAtgE.pdf")
# 
# nums_2_use <- subset(plot.data[,4:107], plot.data$Environment == "lab" & plot.data$Genotype == "AtgE")
# corr_heatmapper(nums_2_use, filename = "cyt_corr_labAtgE.pdf")
# 
# 
# nums_2_use <- subset(plot.data[,4:107], plot.data$Environment == "wild" & plot.data$Genotype == "NOD2")
# corr_heatmapper(nums_2_use, filename = "cyt_corr_wildNOD2.pdf")
# 
# nums_2_use <- subset(plot.data[,4:107], plot.data$Environment == "lab" & plot.data$Genotype == "NOD2")
# corr_heatmapper(nums_2_use, filename = "cyt_corr_labNOD2.pdf")



#
#
#
#### corr by mouse

mouse_heatmapper <- function(input, filename="corr.pdf"){
  
  input <- t(input)
  
  hah2 = corr.test(input, method = "spearman")
  corr_mat = as.matrix(hah2$r)
  
  orig_p = hah2$p
  orig_p[col(orig_p) == row(orig_p) | upper.tri(orig_p)] <- NA
  new_ps = p.adjust(orig_p, method = "BH", n = length(orig_p))
  new_p_mat = matrix(new_ps, ncol = ncol(orig_p), nrow = nrow(orig_p))
  
  new_p_mat[upper.tri(new_p_mat)] = t(new_p_mat)[upper.tri(new_p_mat)]
  diag(new_p_mat) = rep(1,nrow(new_p_mat))
  
  p_list = c(new_p_mat)
  p_dims = dim(new_p_mat)
  
  sym_list = c()
  for (i in p_list){
    if (is.na(i)){
      sym = ""
    }
    else {
      if (i > 0){
        sym = ""
        if (i < 0.05){
          sym = "*"
          if (i<0.01){
            sym = "**"
            if(i<0.001){
              sym = "***"
            }
          }
        }
      }
    }
    sym_list = c(sym_list, sym)
  }
  
  sym_mat = matrix(sym_list, ncol = p_dims[1], nrow = p_dims[2])
  colnames(sym_mat) <- colnames(corr_mat)
  rownames(sym_mat) <- rownames(corr_mat)
  
  corr_mat[is.na(corr_mat)] <- 0
  
  corr_mat2=corr_mat[hclust(dist(corr_mat))$labels[hclust(dist(corr_mat))$order],
                     hclust(dist(corr_mat))$labels[hclust(dist(corr_mat))$order]]
  
  namer2 <- gsub("corr", "dend", filename)
  pdf(namer2, height = 7, width = 10)
  plot(as.dendrogram(hclust(dist(corr_mat))))
  dev.off()
  
  df <- data.frame(names = colnames(corr_mat2))
  df <- merge(plot.data[,1:3], df, by.x = "mouse_id", by.y = "names")
  df$Environment <- gsub("NYU", "lab", df$Location)
  df$Environment <- gsub("SF", "wild", df$Environment)
  df$Environment <- factor(df$Environment, levels=c("lab","wild"))
  df$Genotype <- df$Genetics
  df$Genotype <- gsub("AtgW", "B6", df$Genotype)
  
  #df <- df[order(df$Challenge, df$Cytokine),]
  #corr_mat2 <- corr_mat[as.character(df$og_names),as.character(df$og_names)]
  corr_mat2[upper.tri(corr_mat2, diag = F)] = 0
  
  va <- rowAnnotation(df = df[,4:5],
                      col = list(Environment=c(lab = "mediumorchid3", wild = "red3"),
                                 Genotype=c(AtgE="dodgerblue2", AtgH = "navy",
                                             B6="green4", NOD2="orange")))
  
  ha <- columnAnnotation(df = df[,4:5],
                         col = list(Environment=c(lab = "mediumorchid3", wild = "red3"),
                                    Genotype=c(AtgE="dodgerblue2", AtgH = "navy",
                                               B6="green4", NOD2="orange")))
  
  hah <- Heatmap(corr_mat2, 
                 col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
                 top_annotation = ha,
                 #split = df$Challenge,
                 cluster_columns = F, cluster_rows = F,
                 show_row_names = T, show_column_names = F,
                 row_names_gp = gpar(fontsize = 6),
                 row_names_side = "left",
                 heatmap_legend_param=list(title = "Correlation"))+va
  
  namer = filename
  pdf(namer, height = 8, width = 8)
  draw(hah)
  dev.off()
}

nums_2_use <- plot.data[,4:107]
mouse_heatmapper(nums_2_use, filename = "mouse_corr_all.pdf")

nums_2_use <- subset(plot.data[,4:107], plot.data$Environment == "wild")
mouse_heatmapper(nums_2_use, filename = "mouse_corr_wild.pdf")

nums_2_use <- subset(plot.data[,4:107], plot.data$Environment == "lab")
mouse_heatmapper(nums_2_use, filename = "mouse_corr_lab.pdf")

#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#

##

#
#
#
#
#
#
#



full_cyt <- read.table("full_cytokines.txt", T, '\t')

#Wild

wild_cyt <- subset(full_cyt, Location == "SF")
nums <- data.matrix(wild_cyt[,6:ncol(wild_cyt)])

challenges <- as.character(unique(wild_cyt$Challenge))

thing <- c(1:20)
full_combos <- data.frame(thing)
for (i in 1:length(challenges)){
  c1 <- subset(wild_cyt, Challenge == challenges[i])
  nums_c1 <- data.matrix(c1[,6:ncol(c1)])
  colnames(nums_c1) <- paste(colnames(nums_c1), challenges[i], sep=':')
  full_combos <- cbind(full_combos, nums_c1)
}

nums2 <- full_combos[,-1]
hah2 = corr.test(nums2, method = "spearman")
corr_mat = as.matrix(hah2$r)

orig_p = hah2$p
orig_p[col(orig_p) == row(orig_p) | upper.tri(orig_p)] <- NA
new_ps = p.adjust(orig_p, method = "BH", n = length(orig_p))
new_p_mat = matrix(new_ps, ncol = ncol(orig_p), nrow = nrow(orig_p))

new_p_mat[upper.tri(new_p_mat)] = t(new_p_mat)[upper.tri(new_p_mat)]
diag(new_p_mat) = rep(1,nrow(new_p_mat))

p_list = c(new_p_mat)
p_dims = dim(new_p_mat)

sym_list = c()
for (i in p_list){
  if (is.na(i)){
    sym = ""
  }
  else {
    if (i > 0){
      sym = ""
      if (i < 0.05){
        sym = "*"
        if (i<0.01){
          sym = "**"
          if(i<0.001){
            sym = "***"
          }
        }
      }
    }
  }
  sym_list = c(sym_list, sym)
}

sym_mat = matrix(sym_list, ncol = p_dims[1], nrow = p_dims[2])
colnames(sym_mat) <- colnames(corr_mat)
rownames(sym_mat) <- rownames(corr_mat)


df <- data.frame(names = colnames(corr_mat))
df$Cytokine <- gsub(":.*", "", df$names)
df$Challenge <- gsub(".*:", "", df$names)

va <- columnAnnotation(df = df[,2:3],
                    col = list(Cytokine=c(IL.1a = "red", IL.5 = "blue", IFN.y = "green",
                                          TNF.a = "pink", KC = "purple", MCP.1 = "orange",
                                          IL.1b = "yellow", IL.10 = "turquoise",
                                          CCL3 = "grey", CCL4 = "black", IL.17 = "beige",
                                          IL.6 = "brown", IL.13 = "green4"),
                               Challenge=c(Bacillus="red", Bacteroides = "blue",
                                           Candida="yellow", CD3_CD8.="green",
                                           Clostridium="pink", PBS="purple", Psuedo="orange",
                                           Staph_A="black")))

corr_mat[is.na(corr_mat)] <- 0
hah <- Heatmap(corr_mat, 
               col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
               top_annotation = va,
               cluster_columns = T, cluster_rows = T,
               show_row_names = T, show_column_names = F,
               row_names_gp = gpar(fontsize = 6),
               heatmap_legend_param=list(title = "Correlation"))

pdf("cyt_corr_wild.pdf", height = 8, width = 8)
draw(hah)
dev.off()

#Lab

lab_cyt <- subset(full_cyt, Location == "NYU")
nums <- data.matrix(lab_cyt[,6:ncol(lab_cyt)])

challenges <- as.character(unique(lab_cyt$Challenge))

thing <- c(1:20)
full_combos <- data.frame(thing)
for (i in 1:length(challenges)){
  c1 <- subset(lab_cyt, Challenge == challenges[i])
  nums_c1 <- data.matrix(c1[,6:ncol(c1)])
  colnames(nums_c1) <- paste(colnames(nums_c1), challenges[i], sep=':')
  full_combos <- cbind(full_combos, nums_c1)
}

nums2 <- full_combos[,-1]
hah2 = corr.test(nums2, method = "spearman")
corr_mat = as.matrix(hah2$r)

orig_p = hah2$p
orig_p[col(orig_p) == row(orig_p) | upper.tri(orig_p)] <- NA
new_ps = p.adjust(orig_p, method = "BH", n = length(orig_p))
new_p_mat = matrix(new_ps, ncol = ncol(orig_p), nrow = nrow(orig_p))

new_p_mat[upper.tri(new_p_mat)] = t(new_p_mat)[upper.tri(new_p_mat)]
diag(new_p_mat) = rep(1,nrow(new_p_mat))

p_list = c(new_p_mat)
p_dims = dim(new_p_mat)

sym_list = c()
for (i in p_list){
  if (is.na(i)){
    sym = ""
  }
  else {
    if (i > 0){
      sym = ""
      if (i < 0.05){
        sym = "*"
        if (i<0.01){
          sym = "**"
          if(i<0.001){
            sym = "***"
          }
        }
      }
    }
  }
  sym_list = c(sym_list, sym)
}

sym_mat = matrix(sym_list, ncol = p_dims[1], nrow = p_dims[2])
colnames(sym_mat) <- colnames(corr_mat)
rownames(sym_mat) <- rownames(corr_mat)


df <- data.frame(names = colnames(corr_mat))
df$Cytokine <- gsub(":.*", "", df$names)
df$Challenge <- gsub(".*:", "", df$names)

va <- columnAnnotation(df = df[,2:3],
                       col = list(Cytokine=c(IL.1a = "red", IL.5 = "blue", IFN.y = "green",
                                             TNF.a = "pink", KC = "purple", MCP.1 = "orange",
                                             IL.1b = "yellow", IL.10 = "turquoise",
                                             CCL3 = "grey", CCL4 = "black", IL.17 = "beige",
                                             IL.6 = "brown", IL.13 = "green4"),
                                  Challenge=c(Bacillus="red", Bacteroides = "blue",
                                              Candida="yellow", CD3_CD8.="green",
                                              Clostridium="pink", PBS="purple", Psuedo="orange",
                                              Staph_A="black")))

corr_mat[is.na(corr_mat)] <- 0
hah <- Heatmap(corr_mat, 
               col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
               top_annotation = va,
               cluster_columns = T, cluster_rows = T,
               show_row_names = T, show_column_names = F,
               row_names_gp = gpar(fontsize = 6),
               heatmap_legend_param=list(title = "Correlation"))

pdf("cyt_corr_lab.pdf", height = 8, width = 8)
draw(hah)
dev.off()
