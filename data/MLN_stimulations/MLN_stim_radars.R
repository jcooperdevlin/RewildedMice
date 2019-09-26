###### # # # ##  With fixed MLN data

### Make radars new & improved

setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/data/MLN_stimulations")

library(ComplexHeatmap)
library(ggplot2)
library(reshape)
library(circlize)
library(ggsci)
library(ggrepel)
library(gridExtra)

#save time read it in
mln_cyt_keep = read.table("MLN_stimulation_flat.txt", T, '\t')

meta <- read.table("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/data/metadata/mice_metadata.11.19_mouse_id.txt", 
                   header=T, sep='\t')
rownames(meta) <- as.character(meta$mouse_id)
meta_keep <- meta[as.character(mln_cyt_keep$mouse_id),]

pc <- prcomp(log2(mln_cyt_keep[,-1]+1), scale=T)
huh<-pc$rotation

plotter <- data.frame(PC1=pc$x[,1], PC2=pc$x[,2], mln_cyt_keep, meta_keep)

##
stim_melt <- melt(plotter, measure.vars = colnames(plotter)[4:107], id.vars = colnames(plotter)[c(3,109,110,113)])
stim_melt$challenge <- gsub(".*_", "", stim_melt$variable)
stim_melt$cytokine <- gsub("_.*", "", stim_melt$variable)

ggplot(stim_melt, aes(x=Genotype, y=log2(value+1), color=Genotype)) +
  geom_boxplot() + geom_jitter(width=0.1) +
  theme_bw() +
  facet_wrap(~challenge)

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

# remake some radars

stim_melt <- melt(plotter, measure.vars = colnames(plotter)[4:107], id.vars = colnames(plotter)[c(3,109,110,113)])
stim_melt$challenge <- gsub(".*_", "", stim_melt$variable)
stim_melt$cytokine <- gsub("_.*", "", stim_melt$variable)


challenges <- c("BacillusS","BacteroidesV","CandidaA",
                "CD3.CD28","ClostridiumP","PseudomonasA", 
                "StaphA", "PBS")

orderer <- stim_melt[order(stim_melt$value, decreasing = T),]
orderer <- unique(orderer$cytokine)
orderer <- gsub("\\.", "-", orderer)

glist = list()
for (j in 1:length(challenges)){
  full_df2 <- subset(stim_melt, challenge == challenges[j])
  
  nums2 <- log2(full_df2$value+1)
  
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
  #plotter_slim$labelery <- rep(max(plotter_slim$nums), nrow(plotter_slim))
  plotter_slim$labelery <- rep(10, nrow(plotter_slim))
  plotter_slim$labelerx <- ""
  for(i in 1:nrow(plotter_slim)){
    if(plotter_slim$Environment[i]=="wild"){
      plotter_slim$labelerx[i]<-as.character(plotter_slim$cyto[i])
    }
  }
  
  plotter_slim$nums[plotter_slim$nums>10] <- 10
  
  #filenamer <- paste0("radars/radar_plot_", challenges[j], ".pdf")
  #pdf(filenamer)
  glist[[j]] <- ggplot(plotter_slim, aes(x=cyto, y=nums)) +
    geom_polygon(aes(group = Environment, color = Environment,fill = Environment),
                 alpha=0.4, size = 1, show.legend = TRUE) +
    coord_polar(clip='on')+
    scale_y_continuous(limits = c(0,10)) +
    geom_text(data=plotter_slim, aes(x=cyto, y=labelery, label=labelerx), size = 5) + 
    scale_color_manual(values = c("red3", "mediumpurple1"))+
    scale_fill_manual(values = c("red3", "mediumpurple1"))+
    xlab("Cytokines") + ylab("Log2 pg/mL") +
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


pdf("all_radars.pdf", height = 8, width = 22)
grid.arrange(glist[[1]], glist[[2]], glist[[3]],
             glist[[4]], glist[[5]], glist[[6]],
             glist[[7]], glist[[8]], ncol = 4)
dev.off()



challenges <- c("BacillusS","BacteroidesV","CandidaA",
                "CD3.CD28","ClostridiumP","PseudomonasA", 
                "StaphA", "PBS")
glist = list()
for (j in 1:length(challenges)){
  full_df2 <- subset(stim_melt, challenge == challenges[j])
  
  nums2 <- log2(full_df2$value+1)
  
  stuff <- full_df2$Genotype
  cyto <- full_df2$cytokine
  
  plotter <- data.frame(stuff, cyto, nums2)
  
  plotter_slim <- aggregate(plotter$nums2, list(plotter$stuff, plotter$cyto),mean)
  colnames(plotter_slim) <- c("Genotype", "cyto", 'nums')
  
  plotter_slim$Genotype <- factor(plotter_slim$Genotype, levels = c("B6", "NOD2", "AtgE", "AtgH"))
  plotter_slim$cyto <- gsub("\\.", "-", plotter_slim$cyto)

  plotter_slim=plotter_slim[order(plotter_slim$nums, decreasing = T),]
  #plotter_slim$cyto <- factor(plotter_slim$cyto, levels = rev(unique(plotter_slim$cyto)))
  plotter_slim$cyto <- factor(plotter_slim$cyto, levels = rev(orderer))
  plotter_slim=plotter_slim[order(plotter_slim$cyto, decreasing = T),]
  plotter_slim$labelery <- rep(10, nrow(plotter_slim))
  
  plotter_slim$nums[plotter_slim$nums>10] <- 10
  #filenamer <- paste0("radars/radar_plot_", challenges[j], ".pdf")
  #pdf(filenamer)
  glist[[j]] <- ggplot(plotter_slim, aes(x=cyto, y=nums)) +
    geom_polygon(aes(group = Genotype, color = Genotype,fill = Genotype),
                 alpha=0.4, size = 1, show.legend = TRUE) +
    coord_polar(clip='on')+
    scale_y_continuous(limits = c(0,10)) +
    geom_text(data=plotter_slim, aes(x=cyto, y=labelery, label=cyto), size = 4) + 
    #scale_color_manual(values = c("mediumpurple1", "red3"))+
    #scale_fill_manual(values = c("mediumpurple1", "red3"))+
    xlab("Cytokines") + ylab("Log pg/mL") +
    ggtitle(challenges[j]) +
    theme(
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


pdf("all_radars_genotype.pdf", height = 17, width = 42)
grid.arrange(glist[[1]], glist[[2]], glist[[3]],
             glist[[4]], glist[[5]], glist[[6]],
             glist[[7]], glist[[8]], ncol = 4)
dev.off()



challenges <- c("BacillusS","BacteroidesV","CandidaA",
                "CD3.CD28","ClostridiumP","PseudomonasA", 
                "StaphA", "PBS")
glist = list()
for (j in 1:length(challenges)){
  full_df2 <- subset(stim_melt, challenge == challenges[j])
  
  nums2 <- log2(full_df2$value+1)
  
  stuff <- full_df2$Gender
  cyto <- full_df2$cytokine
  
  plotter <- data.frame(stuff, cyto, nums2)
  
  plotter_slim <- aggregate(plotter$nums2, list(plotter$stuff, plotter$cyto),mean)
  colnames(plotter_slim) <- c("Sex", "cyto", 'nums')
  
  plotter_slim$Sex <- factor(plotter_slim$Sex, levels = c("M", "F"))
  plotter_slim$cyto <- gsub("\\.\\.", "_", plotter_slim$cyto)
  plotter_slim$cyto <- gsub("\\.", "_", plotter_slim$cyto)
  plotter_slim$cyto <- substr(plotter_slim$cyto,1,nchar(plotter_slim$cyto)-1)
  plotter_slim=plotter_slim[order(plotter_slim$nums, decreasing = T),]
  plotter_slim$cyto <- factor(plotter_slim$cyto, levels = rev(unique(plotter_slim$cyto)))
  plotter_slim=plotter_slim[order(plotter_slim$cyto, decreasing = T),]
  plotter_slim$labelery <- rep(max(plotter_slim$nums), nrow(plotter_slim))
  plotter_slim$labelerx <- ""
  for(i in 1:nrow(plotter_slim)){
    if(plotter_slim$Sex[i]=="M"){
      plotter_slim$labelerx[i]<-as.character(plotter_slim$cyto[i])
    }
  }
  
  #filenamer <- paste0("radars/radar_plot_", challenges[j], ".pdf")
  #pdf(filenamer)
  glist[[j]] <- ggplot(plotter_slim, aes(x=cyto, y=nums)) +
    geom_polygon(aes(group = Sex, color = Sex,fill = Sex),
                 alpha=0.4, size = 1, show.legend = TRUE) +
    coord_polar(clip='on')+
    geom_text(data=plotter_slim, aes(x=cyto, y=labelery, label=labelerx), size = 5) + 
    scale_color_manual(values = c("dodgerblue2", "purple"))+
    scale_fill_manual(values = c("dodgerblue2", "purple"))+
    xlab("Cytokines") + ylab("Log pg/mL") +
    ggtitle(challenges[j]) +
    theme(
      #plot.margin=margin(10, 10, 5, 5),
      axis.text.x=element_blank(),
      axis.text.y = element_text(size = 12, color = 'black'),
      axis.title = element_text(size = 15),
      legend.text = element_text(size = 12),
      legend.title = element_text(size=15))
  #print(g)
  #dev.off()
}


pdf("all_radars_gender.pdf", height = 12, width = 32)
grid.arrange(glist[[1]], glist[[2]], glist[[3]],
             glist[[4]], glist[[5]], glist[[6]],
             glist[[7]], glist[[8]], ncol = 4)
dev.off()
