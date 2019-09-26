#### MLN plots Figure 2

####
setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/data/MLN_stimulations")

library(ComplexHeatmap)
library(ggplot2)
library(reshape)
library(circlize)
library(gridExtra)

stims <- read.table("MLN_stimulation_flat.txt", T)

stim_melt <- melt(stims, measure.vars = colnames(stims)[-1], id.vars = colnames(stims)[1])
stim_melt$cytokine <- gsub("_.*", "", stim_melt$variable)
stim_melt$stimulation <- gsub(".*_", "", stim_melt$variable)



#
#
#
#
#
#### figure 2A log fold change of each cytokine as compared to PBS
meta <- read.table("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/data/metadata/mice_metadata.11.19_mouse_id.txt", 
                   header=T, sep='\t')
rownames(meta) <- as.character(meta$mouse_id)
cyt_keep <- intersect(unique(stim_melt$mouse_id), meta$mouse_id)
meta <- meta[cyt_keep,]

stims <- unlist(unique(stim_melt$stimulation))
stims=stims[stims!="PBS"]
cyts <- unlist(unique(stim_melt$cytokine))

fc_df <- data.frame(stim_melt[1,], mult=NA, fc=NA)
for (i in 1:length(stims)){
  for(j in 1:length(cyts)){
    df = subset(stim_melt, stimulation == stims[i] & cytokine == cyts[j])
    pbs = subset(stim_melt, stimulation == "PBS" & cytokine == cyts[j])
    
    rownames(df) <- df$mouse_id
    rownames(pbs) <- pbs$mouse_id
    df[df==0]<-1
    pbs[pbs==0]<-1
    
    pbs <- pbs[rownames(df),]
    
    mult <- ((pbs$value-df$value)/pbs$value)
    df$mult <- -1
    df$mult[mult<0]<- 1
    df$fc <- log2((abs((pbs$value-df$value)/pbs$value))+1)
    df$fc <- df$fc*df$mult
    
    fc_df <- rbind(fc_df, df)
  }
}
fc_df = fc_df[-1,]


fc_mean <- aggregate(data = fc_df, fc ~ cytokine + stimulation, "mean")
fc_heat <- cast(fc_mean,  stimulation ~ cytokine)
rownames(fc_heat) <- fc_heat[,1]
fc_heat = data.frame(fc_heat[,-1])

colnames(fc_heat) <- gsub("\\.", "-", colnames(fc_heat))


labs = colnames(fc_heat)
ll <- columnAnnotation(labels = anno_text(labs, #gp=gpar(fontsize=5),
                                          which = "column", rot=45, 
                                          just = 'right', offset=unit(1.6,"cm")), 
                       height = unit(1.6,"cm"))

pdf("fig2A.pdf", height = 5, width = 9)
Heatmap(fc_heat, bottom_annotation = ll,
        show_column_names = F,
        col = colorRamp2(c(0,5), c("white", "red3")),
        heatmap_legend_param = list(title = "Log2\nFoldChange")
        )
dev.off()

fig_2a_heat <- grid.grabExpr(draw(
  Heatmap(fc_heat, bottom_annotation = ll,
          show_column_names = F,
          col = colorRamp2(c(0,5), c("white", "red3")),
          heatmap_legend_param = list(title = "Log2\nFoldChange")
  )
))


fc_mean <- aggregate(data = fc_df, fc ~ cytokine + stimulation + mouse_id, "mean")
fc_heat <- cast(fc_mean,  mouse_id ~ stimulation + cytokine)
rownames(fc_heat) <- fc_heat[,1]
fc_heat = data.frame(fc_heat[,-1])

meta <- meta[rownames(fc_heat),]

hh <- rowAnnotation(df = data.frame(Environment=meta$Environment, Genetics=meta$Genotype),
                    col = list(Environment = c(lab = "mediumorchid3",
                                               wild = "red3"),
                               Genetics = c(AtgW = "darkorange1",
                                            AtgE = "dodgerblue2",
                                            AtgH = "navyblue",
                                            B6 = "darkorange1",
                                            NOD2 = "mediumspringgreen")))

labs = colnames(fc_heat)
ll <- columnAnnotation(labels = anno_text(labs, gp=gpar(fontsize=5),
                                          which = "column", rot=60, 
                                          just = 'right', offset=unit(1.4,"cm")), 
                       height = unit(1.4,"cm"))

pdf("fig2a_indv.pdf", height = 7, width = 11)
Heatmap(fc_heat,
        bottom_annotation = ll,
        show_row_names = F, show_column_names = F,
        col = colorRamp2(c(-1,0,7), c("blue", "white", "red3")),
        heatmap_legend_param = list(title = "Cytokine")
) + hh
dev.off()



#
#
##heat env geno combo

fc_mean <- aggregate(data = fc_df, fc ~ cytokine + stimulation + mouse_id, "mean")
fc_heat <- cast(fc_mean,  mouse_id ~ stimulation + cytokine)
rownames(fc_heat) <- fc_heat[,1]
fc_heat = data.frame(fc_heat[,-1])

meta <- meta[rownames(fc_heat),]

fc_heat <- cbind(fc_heat, meta$Genotype, meta$Environment)
fc_melt <- melt(fc_heat, measure.vars = colnames(fc_heat)[1:91], id.vars = colnames(fc_heat)[92:93])
fc_melt$combo <- paste0(fc_melt$`meta$Environment`, ":", fc_melt$`meta$Genotype`)

fc_mean <- aggregate(data = fc_melt, value ~ combo + variable, "mean")
fc_heat <- cast(fc_mean,  combo ~ variable)
rownames(fc_heat) <- fc_heat[,1]
fc_heat = data.frame(fc_heat[,-1])

labs = colnames(fc_heat)
ll <- columnAnnotation(labels = anno_text(labs, gp=gpar(fontsize=5),
                                          which = "column", rot=60, 
                                          just = 'right', offset=unit(1.7,"cm")), 
                       height = unit(1.7,"cm"))

pdf("fig2a_combo.pdf", height = 5, width = 11)
Heatmap(fc_heat,
        bottom_annotation = ll,
        show_row_names = T, show_column_names = F,
        col = colorRamp2(c(-1,0,5), c("blue", "white", "red3")),
        heatmap_legend_param = list(title = "Cytokine")
)
dev.off()




#
#
#
#
#
#
# fold change radars?
fc_mean <- aggregate(data = fc_df, fc ~ cytokine + stimulation + mouse_id, "mean")
fc_heat <- cast(fc_mean,  mouse_id ~ stimulation + cytokine)
rownames(fc_heat) <- fc_heat[,1]
fc_heat = data.frame(fc_heat[,-1])

meta <- meta[rownames(fc_heat),]

fc_heat <- cbind(fc_heat, meta$Genotype, meta$Environment, meta$mouse_id)
fc_melt <- melt(fc_heat, measure.vars = colnames(fc_heat)[1:91], id.vars = colnames(fc_heat)[92:94])
fc_melt$cytokine <- gsub(".*_", "", fc_melt$variable)
fc_melt$challenge <- gsub("_.*", "", fc_melt$variable)

colnames(fc_melt)[1:3] <- c("Genotype", "Environment", "mouse_id")

challenges <- c("BacillusS","BacteroidesV","CandidaA",
                "CD3.CD28","ClostridiumP","PseudomonasA", 
                "StaphA")

orderer <- fc_melt[order(fc_melt$value, decreasing = T),]
orderer <- unique(orderer$cytokine)
orderer <- gsub("\\.", "-", orderer)

glist = list()
for (j in 1:length(challenges)){
  full_df2 <- subset(fc_melt, challenge == challenges[j])
  
  nums2 <- full_df2$value
  
  stuff <- full_df2$Environment
  cyto <- full_df2$cytokine
  
  plotter <- data.frame(stuff, cyto, nums2)
  
  plotter_slim <- aggregate(plotter$nums2, list(plotter$stuff, plotter$cyto),mean)
  colnames(plotter_slim) <- c("Environment", "cyto", 'nums')
  
  plotter_slim$Environment <- factor(plotter_slim$Environment, levels = c("wild", "lab"))
  plotter_slim$cyto <- gsub("\\.", "-", plotter_slim$cyto)
  
  plotter_slim=plotter_slim[order(plotter_slim$nums, decreasing = T),]
  plotter_slim$cyto <- factor(plotter_slim$cyto, levels = rev(orderer))
  plotter_slim=plotter_slim[order(plotter_slim$cyto, decreasing = T),]
  plotter_slim$labelery <- rep(max(plotter_slim$nums), nrow(plotter_slim))
  #plotter_slim$labelery <- rep(6, nrow(plotter_slim))
  plotter_slim$labelerx <- ""
  for(i in 1:nrow(plotter_slim)){
    if(plotter_slim$Environment[i]=="wild"){
      plotter_slim$labelerx[i]<-as.character(plotter_slim$cyto[i])
    }
  }
  
  #plotter_slim$nums[plotter_slim$nums>6] <- 6
  
  #filenamer <- paste0("radars/radar_plot_", challenges[j], ".pdf")
  #pdf(filenamer)
  glist[[j]] <- ggplot(plotter_slim, aes(x=cyto, y=nums)) +
    geom_polygon(aes(group = Environment, color = Environment,fill = Environment),
                 alpha=0.4, size = 1, show.legend = TRUE) +
    coord_polar(clip='on')+
    #scale_y_continuous(limits = c(0,6)) +
    geom_text(data=plotter_slim, aes(x=cyto, y=labelery, label=labelerx), size = 5) + 
    scale_color_manual(values = c("red3", "mediumpurple1"))+
    scale_fill_manual(values = c("red3", "mediumpurple1"))+
    xlab("Cytokines") + ylab("Log2 FoldChange") +
    ggtitle(challenges[j]) +
    theme(
      panel.background = element_rect(fill='white'),
      #plot.margin=margin(10, 10, 5, 5),
      panel.grid.major = element_line(color='grey92'),
      axis.text.x=element_blank(),
      axis.text.y = element_text(size = 12, color = 'black'),
      axis.title = element_text(size = 15),
      legend.text = element_text(size = 12),
      legend.title = element_text(size=15))
  #print(g)
  #dev.off()
}


pdf("all_radars_fc_7.10.pdf", height = 8, width = 22)
grid.arrange(glist[[1]], glist[[2]], glist[[3]],
             glist[[4]], glist[[5]], glist[[6]],
             glist[[7]], ncol = 4)
dev.off()

###|##
#
#
#
#
#genotype
challenges <- c("BacillusS","BacteroidesV","CandidaA",
                "CD3.CD28","ClostridiumP","PseudomonasA", 
                "StaphA")
glist = list()
for (j in 1:length(challenges)){
  full_df2 <- subset(fc_melt, challenge == challenges[j])
  
  nums2 <- full_df2$value
  
  stuff1 <- full_df2$Genotype
  stuff2 <- full_df2$Environment
  cyto <- full_df2$cytokine
  
  plotter <- data.frame(stuff1, stuff2, cyto, nums2)
  
  plotter_slim <- aggregate(plotter$nums2, list(plotter$stuff1, plotter$stuff2, plotter$cyto),mean)
  colnames(plotter_slim) <- c("Genotype", "Environment", "cyto", 'nums')
  
  plotter_slim$Genotype <- factor(plotter_slim$Genotype, levels = c("B6", "NOD2", "AtgE", "AtgH"))
  plotter_slim$Environment <- factor(plotter_slim$Environment, levels = c("wild", "lab"))
  plotter_slim$cyto <- gsub("\\.", "-", plotter_slim$cyto)
  
  plotter_slim=plotter_slim[order(plotter_slim$nums, decreasing = T),]
  #plotter_slim$cyto <- factor(plotter_slim$cyto, levels = rev(unique(plotter_slim$cyto)))
  plotter_slim$cyto <- factor(plotter_slim$cyto, levels = rev(orderer))
  plotter_slim=plotter_slim[order(plotter_slim$cyto, decreasing = T),]
  plotter_slim$labelery <- rep(6, nrow(plotter_slim))
  
  plotter_slim$nums[plotter_slim$nums>8] <- 8
  #filenamer <- paste0("radars/radar_plot_", challenges[j], ".pdf")
  #pdf(filenamer)
  glist[[j]] <- ggplot(plotter_slim, aes(x=cyto, y=nums)) +
    geom_polygon(aes(group = Environment, color = Environment,fill = Environment),
                 alpha=0.4, size = 1, show.legend = TRUE) +
    coord_polar(clip='on')+
    scale_y_continuous(limits = c(-1,8)) +
    geom_text(data=plotter_slim, aes(x=cyto, y=labelery, label=cyto), size = 4) + 
    scale_color_manual(values = c("red3", "mediumpurple1"))+
    scale_fill_manual(values = c("red3", "mediumpurple1"))+
    xlab("Cytokines") + ylab("Log2 FoldChange") +
    ggtitle(challenges[j]) +
    theme(
      panel.background = element_rect(fill='white'),
      panel.grid.major = element_line(color="gray92"),
      #plot.margin=margin(10, 10, 5, 5),
      panel.grid.major.x = element_line(colour = "grey95"),
      axis.text.x=element_blank(),
      axis.text.y = element_text(size = 12, color = 'black'),
      axis.title = element_text(size = 15),
      legend.text = element_text(size = 12),
      legend.title = element_text(size=15)) + facet_wrap(~Genotype)
  #print(g)
  #dev.off()
}


pdf("all_radars_genotype_fc.pdf", height = 17, width = 42)
grid.arrange(glist[[1]], glist[[2]], glist[[3]],
             glist[[4]], glist[[5]], glist[[6]],
             glist[[7]], ncol = 4)
dev.off()
candida_g_radar <- glist[[3]]


#
#
#
#
## 
# boxplots with decreasing means compare lab and wild


fc_mean <- aggregate(data = fc_df, fc ~ cytokine + stimulation + mouse_id, "mean")
fc_heat <- cast(fc_mean,  mouse_id ~ stimulation + cytokine)
rownames(fc_heat) <- fc_heat[,1]
fc_heat = data.frame(fc_heat[,-1])

meta <- meta[rownames(fc_heat),]

fc_heat <- cbind(fc_heat, meta$Genotype, meta$Environment, meta$mouse_id)
fc_melt <- melt(fc_heat, measure.vars = colnames(fc_heat)[1:91], id.vars = colnames(fc_heat)[92:94])
fc_melt$cytokine <- gsub(".*_", "", fc_melt$variable)
fc_melt$challenge <- gsub("_.*", "", fc_melt$variable)

colnames(fc_melt)[1:3] <- c("Genotype", "Environment", "mouse_id")

challenges <- c("BacillusS","BacteroidesV","CandidaA",
                "CD3.CD28","ClostridiumP","PseudomonasA", 
                "StaphA")

orderer <- fc_melt[order(fc_melt$value, decreasing = T),]
orderer <- unique(orderer$cytokine)
orderer <- gsub("\\.", "-", orderer)

top_cyt_diff <- data.frame(cyt=NA, stim=NA, lab=NA, wild=NA, diff=NA,pvalue=NA, vars=NA)
for (j in 1:length(challenges)){
  full_df2 <- subset(fc_melt, challenge == challenges[j])
  
  nums2 <- full_df2$value
  
  stuff <- full_df2$Environment
  cyto <- full_df2$cytokine
  
  plotter <- data.frame(stuff, cyto, nums2)
  
  ##pause for stats
  stats = c()
  vars=c()
  for(k in 1:length(unique(plotter$cyto))) {
    statter <- subset(plotter, cyto == unique(plotter$cyto)[k])
    stats=c(stats,wilcox.test(nums2 ~ stuff, data=statter)$p.value)
    
    vars=c(vars,var(statter$nums2))
  }
  
  #
  
  plotter_slim <- aggregate(plotter$nums2, list(plotter$stuff, plotter$cyto),mean)
  colnames(plotter_slim) <- c("Environment", "cyto", 'nums')
  
  plotter_slim$Environment <- factor(plotter_slim$Environment, levels = c("wild", "lab"))
  plotter_slim$cyto <- gsub("\\.", "-", plotter_slim$cyto)
  plotter_slim$stim = challenges[j]
  
  plotter_cast <- cast(plotter_slim, cyto + stim ~ Environment, value = "nums")
  diff = plotter_cast$wild-plotter_cast$lab
  adder <- data.frame(cyt=plotter_cast$cyto, stim=plotter_cast$stim, 
                      lab=plotter_cast$lab, wild=plotter_cast$wild, diff=diff, pvalue=stats,vars=vars)
  top_cyt_diff <- rbind(top_cyt_diff, adder)
  
}
top_cyt_diff <- top_cyt_diff[-1,]
top_cyt_diff<- top_cyt_diff[order(top_cyt_diff$diff, decreasing = T),]
top_cyt_diff$cond <- paste0(top_cyt_diff$stim, ":", top_cyt_diff$cyt)  

#pause make variance supp sigure
top_cyt_diff2 <- top_cyt_diff[order(top_cyt_diff$vars, decreasing =T),]
top_cyt_diff2$cond <- factor(top_cyt_diff2$cond, levels = rev(as.character(top_cyt_diff2$cond)))

top_cyt_agg <- aggregate(data=top_cyt_diff2, vars ~ cyt, "mean")
top_cyt_agg=top_cyt_agg[order(top_cyt_agg$vars, decreasing=T),]
top_cyt_agg$cyt <- factor(top_cyt_agg$cyt, levels = rev(as.character(top_cyt_agg$cyt)))

var_plot=ggplot(top_cyt_agg, aes(cyt, 
                                vars)) +
  geom_col(fill='grey10') + coord_flip() +
  ggtitle("Cytokine Variance") +
  scale_color_manual(values="grey10") +
  theme_bw() +
  ylab("Variance") + xlab("Cytokine") +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12, color="black"),
        legend.position = 'none',
        legend.title = element_text(size=15),
        legend.text = element_text(size=12, color="black")
  )

cyt_heat <- cast(top_cyt_diff, cyt ~ stim, value = "vars")
rownames(cyt_heat) <- cyt_heat[,1]
cyt_heat <- cyt_heat[,-1]
cyt_heat <- as.matrix(data.frame(cyt_heat))

hh=grid.grabExpr(draw(
  Heatmap(cyt_heat,
          col = colorRamp2(c(0, 10), c("white", "navy")),
          cluster_columns = T, cluster_rows = T,
          show_row_names = T, show_column_names = T,
          row_names_gp = gpar(fontsize = 10),
          heatmap_legend_param=list(title = "Variance"))
))

pdf("cytokine_variance.pdf", height = 5, width = 10)
grid.arrange(var_plot, hh, nrow=1)
dev.off()

#

top_cyts_df <- top_cyt_diff[1:10,]
top_cyts_df <- top_cyts_df[order(top_cyts_df$wild, decreasing=T),]
top_cyts <- top_cyts_df$cond

fc_df2 <- fc_df
fc_df2$cytokine <- gsub("\\.", "-", fc_df2$cytokine)
fc_df2$cond <- paste0(fc_df2$stimulation, ":", fc_df2$cytokine)

fc_small <- subset(fc_df2, cond %in% top_cyts)
fc_small$Environment <- rep(meta$Environment, 10)

fc_small$cond <- factor(fc_small$cond, levels = rev(top_cyts))

box_2b <- ggplot(fc_small, aes(cond, fc, color=Environment, fill=Environment)) +
  geom_boxplot(alpha=0.5, outlier.shape=NA) +
  #geom_jitter(width = 0.2) +
  coord_flip() +
  scale_color_manual(values = c("mediumorchid3", "red3")) +
  scale_fill_manual(values = c("mediumorchid3", "red3")) +
  theme_bw() +
  ylab("Log2 Fold Change") + xlab("Stimulation:Cytokine") +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12, color="black"),
        legend.position = 'none',
        legend.title = element_text(size=15),
        legend.text = element_text(size=12, color="black")
        )


fig2_row1 <- arrangeGrob(fig_2a_heat, box_2b, nrow=1, widths = c(2,1))

pdf("fig2ab.pdf", height = 5, width = 12)
grid.arrange(fig2_row1)
dev.off()


#
#
#
#
#
####### fig c and d lower panel of 2
source("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/DIST_EffectSize/var_plotter.R")

stim_cast <- read.table("MLN_stimulation_flat.txt", T, '\t')
meta <- read.table("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/data/metadata/mice_metadata.11.19_mouse_id.txt", T, '\t')

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

effectors=cyt_meta[,c(2,3,4,6,8)]
col_var=cyt_meta$Environment
shape_var=cyt_meta$Genotype

gg_cytokines <- var_plotter(input, effectors, col_var, "Environment",
                            shape_var, "Genotype", "MLN Cytokines")

pdf("fig2cd.pdf", height = 5, width = 15)
grid.arrange(gg_cytokines)
dev.off()




#
#
#
####
#add genotype radar candida only!

pdf("fig2cde.pdf", height = 6, width = 30)
grid.arrange(gg_cytokines, candida_g_radar, nrow=1, widths = c(2,1))
dev.off()

fig2_row2 <- arrangeGrob(gg_cytokines, candida_g_radar, nrow=1, widths = c(2,1))
fig2_row2 <- arrangeGrob(gg_cytokines)
#####fig 2 row 1 and 2

fig2_row12 <- arrangeGrob(fig2_row1, fig2_row2, nrow=2)

#
#
#
#

#
#corr lab v lab and wild v wild

stims <- read.table("MLN_stimulation_flat.txt", T)
meta <- read.table("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/data/metadata/mice_metadata.11.19_mouse_id.txt", 
                   header=T, sep='\t')
rownames(meta) <- as.character(meta$mouse_id)
cyt_keep <- intersect(unique(stims$mouse_id), meta$mouse_id)
meta <- meta[cyt_keep,]

### all corr

corr_mat <- cor(t(stims[,-1]))

hh <- rowAnnotation(df = data.frame(Environment=meta$Environment, Genotype=meta$Genotype),
                    col = list(Environment = c(lab = "mediumorchid3",
                                               wild = "red3"),
                               Genotype = c(AtgW = "darkorange1",
                                            AtgE = "dodgerblue2",
                                            AtgH = "navyblue",
                                            B6 = "darkorange1",
                                            NOD2 = "mediumspringgreen")))

pdf("all_corr.pdf", height = 5, width = 7)
Heatmap(corr_mat,
  show_row_names = F, show_column_names = F,
  heatmap_legend_param = list(title = "Correlation"),
  cluster_rows = T, cluster_columns = T, #row_names_side = 'left',
  #row_names_gp = gpar(fontsize=8),
  #row_names_max_width = unit(10,'cm'),
  col = colorRamp2(c(-1,0,1), c("blue", "white", "red"))) + hh
dev.off()



#
#
#
#lab v lab
lab_meta <- subset(meta, Environment == "lab")
lab_stims <- subset(stims, mouse_id %in% lab_meta$mouse_id)

corr_mat <- cor(t(lab_stims[,-1]))

hh <- rowAnnotation(df = data.frame(Genotype=lab_meta$Genotype),
                    col = list(Genotype = c(AtgW = "darkorange1",
                                            AtgE = "dodgerblue2",
                                            AtgH = "navyblue",
                                            B6 = "darkorange1",
                                            NOD2 = "mediumspringgreen")))

pdf("lab_corr.pdf", height = 5, width = 7)
Heatmap(corr_mat,
        show_row_names = F, show_column_names = F,
        heatmap_legend_param = list(title = "Correlation"),
        cluster_rows = T, cluster_columns = T, #row_names_side = 'left',
        #row_names_gp = gpar(fontsize=8),
        #row_names_max_width = unit(10,'cm'),
        col = colorRamp2(c(-1,0,1), c("blue", "white", "red"))) + hh
dev.off()

lab_corr_heat <- grid.grabExpr(draw(
  Heatmap(corr_mat,
          show_row_names = F, show_column_names = F,
          heatmap_legend_param = list(title = "Correlation"),
          cluster_rows = T, cluster_columns = T, #row_names_side = 'left',
          #row_names_gp = gpar(fontsize=8),
          #row_names_max_width = unit(10,'cm'),
          col = colorRamp2(c(-1,0,1), c("blue", "white", "red"))) + hh
))

#
#
#
#
#
#wild v wild
wild_meta <- subset(meta, Environment == "wild")
wild_stims <- subset(stims, mouse_id %in% wild_meta$mouse_id)

corr_mat <- cor(t(wild_stims[,-1]))

hh <- rowAnnotation(df = data.frame(Genotype=wild_meta$Genotype),
                    col = list(Genotype = c(AtgW = "darkorange1",
                                            AtgE = "dodgerblue2",
                                            AtgH = "navyblue",
                                            B6 = "darkorange1",
                                            NOD2 = "mediumspringgreen")))

pdf("wild_corr.pdf", height = 5, width = 7)
Heatmap(corr_mat,
        show_row_names = F, show_column_names = F,
        heatmap_legend_param = list(title = "Correlation"),
        cluster_rows = T, cluster_columns = T, #row_names_side = 'left',
        #row_names_gp = gpar(fontsize=8),
        #row_names_max_width = unit(10,'cm'),
        col = colorRamp2(c(-1,0,1), c("blue", "white", "red"))) + hh
dev.off()

wild_corr_heat <- grid.grabExpr(draw(
  Heatmap(corr_mat,
          show_row_names = F, show_column_names = F,
          heatmap_legend_param = list(title = "Correlation"),
          cluster_rows = T, cluster_columns = T, #row_names_side = 'left',
          #row_names_gp = gpar(fontsize=8),
          #row_names_max_width = unit(10,'cm'),
          col = colorRamp2(c(-1,0,1), c("blue", "white", "red"))) + hh
))



###fig2 row 3

fig2_row3 <- arrangeGrob(lab_corr_heat, wild_corr_heat, nrow=1)

pdf("fig2.pdf", height = 13, width = 12)
grid.arrange(fig2_row1, fig2_row2, fig2_row3, nrow=3,
             heights = c(0.8,0.9,1))
dev.off()


