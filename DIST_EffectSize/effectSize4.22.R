#### Making Principle Coordinate plots and analyzing effect sizes

#### 4.22

library(ggplot2)
library(circlize)
library(gridExtra)
library(reshape)
library(ape)
library(scales)
library(ggrepel)
library(MDMR)
library(ggsci)
library(ggsignif)

setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data")

type_cols <- c(Diarrhea="red2", Environment="dodgerblue2", FACS = "forestgreen",
               Wedge_cage="navyblue",Flow.date="grey50",
               Gender="hotpink",Genotype="darkorange1", Pregnant="purple", WeightGain="forestgreen")


#function for plotting PCoA with Effect Sizes
var_plotter <- function(input, effectors, col_var, col_var_name="Environment", col_var_col=NA,
                        shape_var=NA, shape_var_name=NA, Feature, dist_panel=F){
  
  D <- dist(input, method = "euclidean")
  res <- pcoa(D)
  plotter <- data.frame(res$vectors, col_var=col_var, 
                        shape_var = shape_var)
  
  #good_plotter <- subset(plotter, Axis.1 > quantile(plotter$Axis.1, 0.05) &
  #                         Axis.1 < quantile(plotter$Axis.1, 0.95))
  #good_plotter <- subset(good_plotter, Axis.2 > quantile(plotter$Axis.2, 0.05) &
  #                         Axis.2 < quantile(plotter$Axis.2, 0.95))

  if(is.na(col_var_col)){col_var_col = c(lab="mediumpurple1", wild="red3")}
  
  g0=ggplot(plotter, aes(Axis.1, Axis.2, color=col_var, shape=shape_var)) + 
    geom_point(size=5) +
    guides(color=guide_legend(title=col_var_name), 
           shape=guide_legend(title=shape_var_name)) +
    ylab(paste0("PCo2 ", round(res$values$Broken_stick[2]*100,2), "% expl. variation")) +
    xlab(paste0("PCo1 ", round(res$values$Broken_stick[1]*100,2), "% expl. variation")) +
    scale_color_manual(values = col_var_col) +
    ggtitle(Feature)+
    theme_bw() +
    theme(axis.text = element_text(size=12, color="black"),
          axis.title = element_text(size=13, color="black"),
          legend.text = element_text(size=12, color="black"),
          legend.title = element_text(size=13, color="black"),
          plot.title = element_text(size=14, color="black", 
                                    face='bold',hjust = 0.5))
  
  pr.coo=res$vectors
  plot.axes=c(1,2)
  n <- nrow(input)
  points.stand <- scale(pr.coo[,plot.axes])
  S <- cov(input, points.stand)
  U <- S %*% diag((res$values$Eigenvalues[plot.axes]/(n-1))^(-0.5))
  colnames(U) <- colnames(pr.coo[,plot.axes])
  
  PC = res
  data <- data.frame(obsnames=row.names(PC$vectors), PC$vectors[,1:2])
  datapc <- data.frame(varnames=rownames(U), U*100)
  datapc$var1 <- rescale(datapc$Axis.1, c(min(data$Axis.1),max(data$Axis.1)))
  datapc$var2 <- rescale(datapc$Axis.1, c(min(data$Axis.2),max(data$Axis.2)))
  
  datapc$mult <- abs(datapc$Axis.1*datapc$Axis.2)
  datapc <- datapc[order(datapc$mult, decreasing = T),]
  datapc2 = datapc
  datapc2 = datapc[1:12,]
  
  g_seg1=ggplot(plotter, aes(Axis.1, Axis.2, color=col_var, shape=shape_var)) +
    geom_point(size=5) +
    guides(color=guide_legend(title=col_var_name), 
           shape=guide_legend(title=shape_var_name)) +
    ylab(paste0("PCo2 ", round(res$values$Broken_stick[2]*100,2), "% expl. variation")) +
    xlab(paste0("PCo1 ", round(res$values$Broken_stick[1]*100,2), "% expl. variation")) +
    scale_color_manual(values = col_var_col) +
    ggtitle(Feature)+
    theme_bw() +
    theme(axis.text = element_text(size=12, color="black"),
          axis.title = element_text(size=13, color="black"),
          legend.text = element_text(size=12, color="black"),
          legend.title = element_text(size=13, color="black"),
          plot.title = element_text(size=14, color="black", 
                                    face='bold',hjust = 0.5))
  g_seg=ggplot(plotter, aes(Axis.1, Axis.2)) +
    geom_point(size=3, alpha=0) + #ggtitle("Loadings for Gated Populations") +
    scale_color_manual(values = col_var_col) + coord_equal() + 
    geom_segment(data=datapc2, aes(x=0, y=0, xend=Axis.1, yend=Axis.2), 
                 arrow=arrow(length=unit(0.2,"cm")), alpha=0.5) +
    geom_label_repel(data=datapc2, aes(x=Axis.1, y=Axis.2, label=varnames), 
                     size = 4, force=4, segment.alpha=0.5) +
    theme_void()
  #g_seg
  
  
  
  pc1_hist <- ggplot(plotter, aes(Axis.1, color = col_var, fill= col_var)) +
    geom_density(alpha=0.6) + 
    scale_y_reverse() +
    scale_color_manual(values = col_var_col) +
    scale_fill_manual(values = col_var_col) +
    theme_void()
  
  pc2_hist <- ggplot(plotter, aes(Axis.2, color = col_var, fill= col_var)) +
    geom_density(alpha=0.6) + 
    coord_flip() +
    scale_y_reverse() +
    scale_color_manual(values = col_var_col) +
    scale_fill_manual(values = col_var_col) +
    theme_void()
  
  pc_fake <- ggplot(plotter, aes(Axis.2, color = col_var, fill= col_var)) +
    geom_density(alpha=0.6) + 
    coord_flip() +
    scale_y_reverse() +
    scale_color_manual(values = c("white", "white", "white")) +
    scale_fill_manual(values = c("white", "white", "white")) +
    theme_void()
  
  top_g=arrangeGrob(
    pc2_hist+theme(legend.position = 'none'),
    g0+theme_void()+
      theme(legend.position = 'none',
            plot.title = element_text(size=14, color="black", 
                                      face='bold',hjust = 0.5)), nrow=1, widths = c(1,3))
  
  g_seg_plot <- arrangeGrob(g_seg,pc_fake+theme(legend.position = 'none'),nrow=2, heights = c(3,1))
  #g_seg_plot <- arrangeGrob(pc_fake+theme(legend.position = 'none'),g_seg_plot,nrow=1, widths = c(1,7))
  
  bottom_g = arrangeGrob(
    pc_fake+theme(legend.position = 'none'),
    pc1_hist+theme(legend.position = 'none'), nrow = 1, widths = c(1,3))
  
  pcoa_plot=arrangeGrob(top_g, bottom_g, heights = c(3,1))

  
  
  mdmr.res <- mdmr(X = effectors, D = D)
  es_df=mdmr.res$stat
  es_df$Variable=gsub("1", "", rownames(es_df))
  es_df$Variable=gsub("Weight_gain\\.", "WeightGain", es_df$Variable)
  es_df <- es_df[-1,]
  
  #delta_res=delta(effectors, Y = inputs, dtype = "euclidean", niter = 1, seed = 12345, plot.res = F)
  #delta_res <- data.frame(t(delta_res))
  
  es_df <- es_df[order(es_df$stat, decreasing = T),]
  es_df$Variable <- factor(es_df$Variable, levels = es_df$Variable)
  es_df <- subset(es_df, !is.na(stat))
  es_df <- es_df[order(abs(es_df$stat), decreasing = T),]
  es_df$Variable <- factor(es_df$Variable, levels = es_df$Variable)
  gg_combo=ggplot(es_df, aes(Variable, abs(stat), fill = Variable)) +
    geom_col() +
    scale_x_discrete(limits = rev(levels(es_df$Variable)))+
    guides(fill=guide_legend(title="Feature"))+
    ylab("EffectSize")+xlab("")+
    scale_fill_manual(values=type_cols)+
    theme_bw()+
    coord_flip()+
    theme(legend.position = 'none',
          axis.text = element_text(size=12, color="black"),
          #axis.text = element_blank(),
          #axis.title = element_text(size=13, color="black"),
          axis.title = element_blank(),
          legend.text = element_text(size=12, color="black"),
          legend.title = element_text(size=13, color="black"),
          plot.title = element_text(size=14, face='bold'))
  
  if(dist_panel==T){
    d_mat <- as.matrix(D)
    d_df <- data.frame(
      id1 = rep(colnames(d_mat), ncol(d_mat)),
      id2 = rep(rownames(d_mat), each=ncol(d_mat)),
      var1 = rep(col_var, ncol(d_mat)),
      var2 = rep(col_var, each=ncol(d_mat)),
      value = c(d_mat)
    )
    d_df$lab <- paste0(d_df$var1, ":", d_df$var2)
    d_df$lab <- gsub("wild:lab", "lab:wild", d_df$lab)
    
    d_df <- subset(d_df, value > 0)
    d_df$lab2 <- factor(d_df$lab, levels = c("lab:lab", "wild:wild", "lab:wild"))
    dist_plot=ggplot(d_df, aes(lab2, value, color=lab, fill=lab)) +
      #geom_jitter(width = 0.3, alpha=0.01) +
      #geom_violin(alpha=0.5, outlier.shape = NA) +
      geom_boxplot(alpha=0.2) +
      scale_fill_manual(values = c("mediumorchid3", "red3", "navy")) +
      scale_color_manual(values = c("mediumorchid3", "red3", "navy")) +
      xlab("Comparison") + ylab("Distance") +
      theme_bw() + geom_signif(test = "wilcox.test", 
                               comparisons = combn(unique(d_df$lab),2, simplify = F),
                               y_position = c(max(d_df$value)+0.5, max(d_df$value)+1, max(d_df$value)+1.5),
                               color='black', map_signif_level = T) +
      theme(legend.position='none',
            axis.title = element_text(size=15),
            axis.text = element_text(size=12, color='black'))
    
    return(arrangeGrob(pcoa_plot, g_seg_plot,dist_plot, gg_combo, ncol=4, widths=c(2,2,1.5,1.5)))
    
  } else {
    return(arrangeGrob(pcoa_plot, g_seg_plot, gg_combo, ncol=3, widths=c(2,2,1)))
  }
  
}

## FACS Blood

facs_blood <- read.table("int/data/FACS_data/BLOOD_lymph_FACS_metadata_names.txt", T, '\t')
name_change <- read.table("int/data/FACS_data/lymph_name_change.txt")
colnames(facs_blood)[13:27] <- as.character(name_change$V2)
rownames(facs_blood) <- as.character(facs_blood$mouse_id)

input=log2(facs_blood[,13:27]+1)
rownames(input) <- as.character(facs_blood$mouse_id)

rowZeros <- apply(input, 1, function(x) {sum(x == 0)})
colZeros <- apply(input, 2, function(x) {sum(x == 0)})

input = input[which(rowZeros<0.5*ncol(input)),which(colZeros<0.5*nrow(input)) ]

facs_blood <- facs_blood[rownames(input),]
effectors=facs_blood[,c(4,5,6,8,10)]
col_var=facs_blood$Environment
shape_var=facs_blood$Genotype

gg_blood_lymph <- var_plotter(input, effectors, col_var, "Environment", col_var_col = NA,
                  shape_var, "Genotype", "Blood Lymph Panel", dist_panel=T)

#grid.arrange(gg_blood_lymph)
### fig 1d?
png("int/DIST_EffectSize/Fig_1d.png", 
    height = 5, width = 20, units = 'in', res=300)
grid.arrange(gg_blood_lymph)
dev.off()

pdf("int/DIST_EffectSize/Fig_1d.pdf", 
    height = 5, width = 20)
grid.arrange(gg_blood_lymph)
dev.off()



## FACS MLN

facs_mln <- read.table("int/data/FACS_data/MLN_lymph_FACS_metadata_names.txt", T, '\t')
name_change <- read.table("int/data/FACS_data/lymph_name_change.txt")
colnames(facs_mln)[13:27] <- as.character(name_change$V2)
rownames(facs_mln) <- as.character(facs_mln$mouse_id)

facs_mln <- facs_mln[!is.na(facs_mln$Weight_gain.),]
facs_mln <- facs_mln[!is.na(facs_mln$Total_T_cells),]

input=log2(facs_mln[,13:27]+1)
rownames(input) <- as.character(facs_mln$mouse_id)

rowZeros <- apply(input, 1, function(x) {sum(x == 0)})
colZeros <- apply(input, 2, function(x) {sum(x == 0)})

input = input[which(rowZeros<0.5*ncol(input)),which(colZeros<0.5*nrow(input)) ]

facs_mln <- facs_mln[rownames(input),]

effectors=facs_mln[,c(4,5,6,8,10)]
col_var=facs_mln$Environment
shape_var=facs_mln$Genotype

gg_mln_lymph <- var_plotter(input, effectors, col_var, "Environment", col_var_col=NA,
                  shape_var, "Genotype", "MLN Lymph Panel")

col_var=facs_mln$Flow.date
gg_mln_lymph_date <- var_plotter(input, effectors, col_var, "Week", 
                                 col_var_col=c("dodgerblue2", "green4", "orange"),
                                 shape_var, "Genotype", "MLN Lymph Panel")

pdf("int/DIST_EffectSize/MLN_lymph_pcoa_ES_combo+date.pdf", height = 10, width = 15)
grid.arrange(gg_mln_lymph, gg_mln_lymph_date, nrow=2)
dev.off()


### mln lymph lab only
lab_mice <- as.character(subset(facs_mln, Environment == "lab")$mouse_id)
input2 <- input[lab_mice,]

facs_mln2 <- facs_mln[lab_mice,]

effectors=facs_mln2[,c(7,6)]
effectors$Wege_cage <- factor(effectors$Wege_cage)
colnames(effectors)[1] <- "Wedge_cage"
col_var=facs_mln2$Environment
shape_var=facs_mln2$Genotype

gg_mln_lymph_lab <- var_plotter(input2, effectors, col_var, "Environment",
                            shape_var, "Genotype", "MLN Lymph Panel")
grid.arrange(gg_mln_lymph_lab)


##
#
#### mln lymph wild only
wild_mice <- as.character(subset(facs_mln, Environment == "wild")$mouse_id)
input2 <- input[wild_mice,]

facs_mln2 <- facs_mln[wild_mice,]

effectors <- data.frame(Flow.date=rep(1,101))
effectors$Flow.date[facs_mln2$Flow.date=="10-Aug"] <- 2

effectors=cbind(effectors, facs_mln2[,c(4,6,7,8,9,10)])
#effectors$Wege_cage <- factor(effectors$Wege_cage)
colnames(effectors)[4] <- "Wedge_cage"
col_var=facs_mln2$Environment
shape_var=facs_mln2$Genotype

gg_mln_lymph_wild <- var_plotter(input2, effectors, col_var, "Environment",
                                shape_var, "Genotype", "MLN Lymph Panel")
grid.arrange(gg_mln_lymph_wild)

pdf("int/DIST_EffectSize/MLN_lymph_pcoa_ES_combo+wild.pdf", height = 10, width = 15)
grid.arrange(gg_mln_lymph, gg_mln_lymph_wild, nrow=2)
dev.off()

#
#
#
#### supp figure (investigate separation)
## FACS MLN My

facs_mln <- read.table("int/data/FACS_data/MLN_myeloid_FACS_metadata_names.txt", T, '\t')
rownames(facs_mln) <- as.character(facs_mln$mouse_id)

facs_mln <- facs_mln[!is.na(facs_mln$Weight_gain.),]
facs_mln <- facs_mln[!is.na(facs_mln$Neutrophils),]

input=log2(facs_mln[,13:18]+1)
rownames(input) <- as.character(facs_mln$mouse_id)

rowZeros <- apply(input, 1, function(x) {sum(x == 0)})
colZeros <- apply(input, 2, function(x) {sum(x == 0)})

input = input[which(rowZeros<0.5*ncol(input)),which(colZeros<0.5*nrow(input)) ]

facs_mln <- facs_mln[rownames(input),]

effectors=facs_mln[,c(4,5,6,8,9,10)]
col_var=facs_mln$Environment
shape_var=facs_mln$Genotype

gg_mln_myeloid <- var_plotter(input, effectors, col_var, "Environment",
                            shape_var, "Genotype", "MLN Myeloid Panel")


#pdf('figures/effectSizes.pdf', height = 18, width = 20)
#grid.arrange(gg_blood_lymph,gg_mln_lymph,gg_mln_myeloid,nrow=3)
#dev.off()




#
#

## combo

blood_lymph <- read.table("int/data/FACS_data/BLOOD_lymph_FACS_metadata_names.txt", header=T, sep='\t')
name_change <- read.table("int/data/FACS_data/lymph_name_change.txt")
colnames(blood_lymph)[13:27] <- as.character(name_change$V2)
colnames(blood_lymph)[13:27] <- paste0("BLOOD_lymph_", colnames(blood_lymph[,13:27]))
blood_lymph <- blood_lymph[order(blood_lymph$mouse_id),]
mln_lymph <- read.table("int/data/FACS_data/MLN_lymph_FACS_metadata_names.txt", header=T, sep='\t')
name_change <- read.table("int/data/FACS_data/lymph_name_change.txt")
colnames(mln_lymph)[13:27] <- as.character(name_change$V2)
colnames(mln_lymph)[13:27] <- paste0("MLN_lymph_", colnames(mln_lymph[,13:27]))
mln_lymph <- mln_lymph[order(mln_lymph$mouse_id),]
mln_myeloid <- read.table("int/data/FACS_data/MLN_myeloid_FACS_metadata_names.txt", header=T, sep='\t')
colnames(mln_myeloid)[13:18] <- paste0("MLN_myeloid_", colnames(mln_myeloid[,13:18]))
mln_myeloid <- mln_myeloid[order(mln_myeloid$mouse_id),]


full_facs <- cbind(blood_lymph, mln_lymph[,13:27], mln_myeloid[,13:18])
num_check <- which(is.na(rowSums(full_facs[,13:48])))
full_facs <- full_facs[-num_check,]
rownames(full_facs) <- as.character(full_facs$mouse_id)

input=log2(full_facs[,13:48]+1)
rownames(input) <- as.character(full_facs$mouse_id)

rowZeros <- apply(input, 1, function(x) {sum(x == 0)})
colZeros <- apply(input, 2, function(x) {sum(x == 0)})

input = input[which(rowZeros<0.5*ncol(input)),which(colZeros<0.5*nrow(input)) ]

full_facs <- full_facs[rownames(input),]

effectors=full_facs[,c(4,5,6,8,9,10)]
col_var=full_facs$Environment
shape_var=full_facs$Genotype

gg_facs_all <- var_plotter(input, effectors, col_var, "Environment",
                              shape_var, "Genotype", "All FACS")


pdf('figures/effectSizes.pdf', height = 24, width = 20)
grid.arrange(gg_blood_lymph,gg_mln_lymph,gg_mln_myeloid,gg_facs_all,nrow=4)
dev.off()


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

## MLN stimulation Cytokines
stim_cast <- read.table("int/data/MLN_stimulations/MLN_stimulation_flat.txt", T, '\t')
meta <- read.table("int/data/metadata/mice_metadata.11.19_mouse_id.txt", T, '\t')

nums <- stim_cast[,2:ncol(stim_cast)]
cc <- colnames(nums)
rr <- as.character(stim_cast$mouse_id)
#nums <- scale(nums)
nums <- data.frame(data.matrix(nums))
colnames(nums) <- cc
rownames(nums) <- as.character(stim_cast$mouse_id)

rownames(meta) <- as.character(meta$mouse_id)
cyt_keep = intersect(meta$mouse_id, stim_cast$mouse_id)
cyt_meta = meta[cyt_keep,]

nums <- nums[cyt_keep,]

input=log2(nums+1)

rowZeros <- apply(input, 1, function(x) {sum(x == 0)})
colZeros <- apply(input, 2, function(x) {sum(x == 0)})

input = input[which(rowZeros<0.5*ncol(input)),which(colZeros<0.5*nrow(input)) ]

cyt_meta<-cyt_meta[rownames(input),]

effectors=cyt_meta[,c(2,3,4,6,7,8)]
col_var=cyt_meta$Environment
shape_var=cyt_meta$Genotype

gg_cytokines <- var_plotter(input, effectors, col_var, "Environment", col_var_col=NA,
                           shape_var, "Genotype", "MLN Cytokines")

#pdf("mln_cytokine_pcoa_ES.pdf", height = 5, width = 15)
#grid.arrange(gg_ser_cyt)
#dev.off()


## quick check
checker = input
checker$mouse_id <- rownames(checker)
checker <- data.frame(checker, cyt_meta)
gg_melt <- melt(checker, measure.vars = c("IL.10_ClostridiumP", "IL.10_CD3.CD28",
                                        "IL.10_StaphA", "IL.10_CandidaA", "IFN.y_StaphA",
                                        "IFN.y_PseudomonasA", "IFN.y_StaphA"), 
                id.vars = c("mouse_id", "Genotype", "Environment", "Pregnant"))

g1=ggplot(gg_melt, aes(Genotype, value, color=Genotype)) +
  geom_boxplot() + geom_jitter(width = 0.1)

g1+facet_wrap(~variable, nrow=3)

g2=ggplot(gg_melt, aes(Pregnant, value, color=Pregnant)) +
  geom_boxplot() + geom_jitter(width = 0.1)

g2+facet_wrap(~variable, nrow=3)



colors_clusters = c(pal_d3("category10")(10), pal_d3("category20b")(20), pal_igv("default")(51))
g3=ggplot(gg_melt, aes(Environment:Genotype, value, color=Environment:Genotype)) +
  geom_boxplot() + geom_jitter(width = 0.1) +
  scale_color_manual(values=colors_clusters)

g3+facet_wrap(~variable, nrow=3)

pdf("figures/top_comp.pdf", height=15,width=22)
grid.arrange(
  arrangeGrob(g1+facet_wrap(~variable, nrow=3),
             g2+facet_wrap(~variable, nrow=3),
             ncol=2),
  g3+facet_wrap(~variable, nrow=2), nrow=2)
dev.off()


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
### serum cytokines
ser_cyt <- read.table("int/data/Serum_stimulations/Serum_names_final.txt", T, '\t')
meta <- read.table("int/data/metadata/mice_metadata.11.19_mouse_id.txt", T, '\t')

cyt_keep <- intersect(ser_cyt$sample, meta$mouse_id)

rownames(ser_cyt) <- ser_cyt$sample
rownames(meta) <- meta$mouse_id

ser_cyt_keep <- ser_cyt[cyt_keep,]
meta_keep <- meta[cyt_keep,]

input=log2(ser_cyt_keep[,-1]+1)

rowZeros <- apply(input, 1, function(x) {sum(x == 0)})
colZeros <- apply(input, 2, function(x) {sum(x == 0)})

input = input[which(rowZeros<0.5*ncol(input)),which(colZeros<0.5*nrow(input)) ]

meta_keep<-meta_keep[rownames(input),]

effectors=meta_keep[,c(2,3,4,6,7,8)]
col_var=meta_keep$Environment
shape_var=meta_keep$Genotype

gg_ser_cyt <- var_plotter(input, effectors, col_var, "Environment",
                            shape_var, "Genotype", "Serum Cytokines")

grid.arrange(gg_ser_cyt)

#

#
#
pdf('figures/effectSizes2.pdf', height = 18, width = 35)
grid.arrange(gg_blood_lymph,gg_mln_lymph,
             gg_mln_myeloid,gg_facs_all,
             gg_cytokines, gg_ser_cyt,
             nrow=3, ncol=2)
dev.off()