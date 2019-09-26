### Plot figure 1b and c indiv umaps and some boxplots
### 6.25

setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/FACS_umaps")

library(ggplot2)
library(reshape)
library(ggsci)
library(gridExtra)

colors_clusters = c(pal_d3("category10")(10), pal_d3("category20b")(20), pal_igv("default")(51))
ramper_basic = colorRampPalette(c("grey99","purple4"))(2)
ramper_more = colorRampPalette(c("white","darkorange", "red2", "purple4"))(100)

plotter2D <- function(input_df, type, x, y, num_cols, xlims, ylims, splitter, ramper=NA) {
  plot_plist <- list()
  counter=1
  for(i in num_cols){
    color_item <- input_df[,i]
    color_item <- log2(color_item+1)
    Log2 <- color_item
    
    lab_df <- subset(input_df, Environment == "lab")
    new_lab <- paste0("Lab ", round(length(which(lab_df[,i]>0))/nrow(lab_df)*100,3), "%")
    wild_df <- subset(input_df, Environment == "wild")
    new_wild <- paste0("Wild ", round(length(which(wild_df[,i]>0))/nrow(wild_df)*100,3), "%")
    
    input_df$Environment2<-input_df$Environment
    input_df$Environment2 <- gsub("lab", new_lab, input_df$Environment2)
    input_df$Environment2 <- gsub("wild", new_wild, input_df$Environment2)
    
    if(is.na(ramper)){
      ramper = colorRampPalette(c("grey90", "purple4"))(2)
    } else {ramper=ramper}
    g=ggplot(input_df, aes(input_df[,x], input_df[,y], color = Log2))+
      geom_point(size=0.001, aes(alpha=Log2))+
      ylim(ylims)+xlim(xlims)+
      #ylab(paste0(type, "_2")) + xlab(paste0(type, "_1"))+
      scale_colour_gradientn(colors = ramper)+
      #scale_colour_gradient2(low="gray85", mid="red2", high="purple4")+
      ggtitle(paste0(colnames(input_df[i]))) +
      theme_void() + 
      theme(legend.position='none',
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            plot.title = element_text(hjust = 0.5))
      
    if(splitter==T){
      plot_plist[[counter]]<-ggplotGrob(g+facet_wrap(~Environment2, nrow=2))
    } else {
      plot_plist[[counter]]<-ggplotGrob(g)
    }
    counter=counter+1
  }
  return(plot_plist)
}



#
#
#

#
#
#
#
#
#Blood
# Row 1. TSNE vs UMAP & CD3, CD19, CD4, CD8 Combined

### Blood all

umap_blood_df <- read.table("inputs/umap_combo_Blood.txt", F, '\t')
meta <- read.table("mice_metadata.11.19_mouse_id.txt", T, '\t')
orig_df <- read.table("inputs/Blood_df.txt", F, '\t')
ids <- orig_df$V1

names <- colnames(read.table("name_change.csv", T, ","))

colnames(orig_df) <- c("id", names)

### add metadata

orig_df$id <- factor(orig_df$id, levels = unique(orig_df$id))
orig_df <- orig_df[order(orig_df$id),]

rownames(meta) <- meta$mouse_id
meta <- meta[levels(orig_df$id),]

uniq_ids <- unique(orig_df$id)
orig_df2 <- data.frame(orig_df[1,], Genotype = NA, Environment = NA,
                       Wedge_cage = NA, Gender = NA, Pregnant = NA, Diarrhea = NA, Flow.date=NA)
for (j in 1:length(uniq_ids)){
  curr <- subset(orig_df, id == uniq_ids[j])
  meta_curr <- subset(meta, mouse_id == as.character(uniq_ids[j]))
  
  curr$Genotype <- rep(meta_curr$Genotype, each = nrow(curr))
  curr$Environment <- rep(meta_curr$Environment, each =nrow(curr))
  curr$Wedge_cage <- rep(meta_curr$Wedge_cage, each =nrow(curr))
  curr$Gender <- rep(meta_curr$Gender, each =nrow(curr))
  curr$Pregnant <- rep(meta_curr$Pregnant, each =nrow(curr))
  curr$Diarrhea <- rep(meta_curr$Diarrhea, each =nrow(curr))
  curr$Flow.date <- rep(meta_curr$Flow.date, each =nrow(curr))
  
  orig_df2 <- rbind(orig_df2, curr)
}
orig_df <- orig_df2[-1,]

colnames(orig_df)

# blood all

orig_df$umap_1 <- umap_blood_df$V1
orig_df$umap_2 <- umap_blood_df$V2

umap_xlims <- c(floor(min(orig_df$umap_1)), ceiling(max(orig_df$umap_1)))
umap_ylims <- c(floor(min(orig_df$umap_2)), ceiling(max(orig_df$umap_2)))


thresh <- read.table("Distributions/thresholds/Blood_major_thresholds.txt", T, "\t")
orig_df$CD3[orig_df$CD3<thresh$CD3[1]] <- 0
orig_df$CD19[orig_df$CD19<thresh$CD19[1]] <- 0
orig_df$CD4[orig_df$CD4<thresh$CD4[1]] <- 0
orig_df$CD8[orig_df$CD8<thresh$CD8[1]] <- 0

int_cols <- c(which(colnames(orig_df)=="CD19"))

cd19_blood <- plotter2D(orig_df, "umap", "umap_1", "umap_2", 
                              int_cols, umap_xlims, umap_ylims, splitter=F)

cd19 <- arrangeGrob(grobs=cd19_blood, nrow=1)

int_cols <- c(which(colnames(orig_df)=="CD4"))

cd4_blood <- plotter2D(orig_df, "umap", "umap_1", "umap_2", 
                        int_cols, umap_xlims, umap_ylims, splitter=F)

cd4 <- arrangeGrob(grobs=cd4_blood, nrow=1)



#
#
#
#
#

#
#
#
### Blood CD19-CD44
### Blood CD4-CD62L

reader19 <- paste0("inputs/CD19_umap_combo_Blood.txt")
umap_blood_df <- read.table(reader19, F, '\t')
meta <- read.table("mice_metadata.11.19_mouse_id.txt", T, '\t')
reader19 <- paste0("inputs/CD19_Blood_df.txt")
orig_df <- read.table(reader19, F, '\t')
ids <- orig_df$V1

names <- colnames(read.table("name_change.csv", T, ","))

colnames(orig_df) <- c("id", names)

### add metadata

orig_df$id <- factor(orig_df$id, levels = unique(orig_df$id))
orig_df <- orig_df[order(orig_df$id),]

rownames(meta) <- meta$mouse_id
meta <- meta[levels(orig_df$id),]

uniq_ids <- unique(orig_df$id)
orig_df2 <- data.frame(orig_df[1,], Genotype = NA, Environment = NA,
                       Wedge_cage = NA, Gender = NA, Pregnant = NA, Diarrhea = NA, Flow.date=NA)
for (j in 1:length(uniq_ids)){
  curr <- subset(orig_df, id == uniq_ids[j])
  meta_curr <- subset(meta, mouse_id == as.character(uniq_ids[j]))
  
  curr$Genotype <- rep(meta_curr$Genotype, each = nrow(curr))
  curr$Environment <- rep(meta_curr$Environment, each =nrow(curr))
  curr$Wedge_cage <- rep(meta_curr$Wedge_cage, each =nrow(curr))
  curr$Gender <- rep(meta_curr$Gender, each =nrow(curr))
  curr$Pregnant <- rep(meta_curr$Pregnant, each =nrow(curr))
  curr$Diarrhea <- rep(meta_curr$Diarrhea, each =nrow(curr))
  curr$Flow.date <- rep(meta_curr$Flow.date, each =nrow(curr))
  
  orig_df2 <- rbind(orig_df2, curr)
}
orig_df <- orig_df2[-1,]

colnames(orig_df)

# blood all

orig_df$umap_1 <- umap_blood_df$V1
orig_df$umap_2 <- umap_blood_df$V2

umap_xlims <- c(floor(min(orig_df$umap_1)), ceiling(max(orig_df$umap_1)))
umap_ylims <- c(floor(min(orig_df$umap_2)), ceiling(max(orig_df$umap_2)))

reader <- paste0("Distributions/thresholds/CD19_Blood_minor_thresholds.txt")
thresh <- read.table(reader, T, "\t")

orig_df$CD43_1B11[orig_df$CD43_1B11<thresh$CD43_1B11[1]] <- 0
orig_df$PD1[orig_df$PD1<thresh$PD1[1]] <- 0
orig_df$CD25[orig_df$CD25<thresh$CD25[1]] <- 0
orig_df$CD44[orig_df$CD44<thresh$CD44[1]] <- 0
orig_df$CD127[orig_df$CD127<thresh$CD127[1]] <- 0
orig_df$CXCR3[orig_df$CXCR3<thresh$CXCR3[1]] <- 0
orig_df$KLRG1[orig_df$KLRG1<thresh$KLRG1[1]] <- 0
orig_df$CD27[orig_df$CD27<thresh$CD27[1]] <- 0
orig_df$CD69[orig_df$CD69<thresh$CD69[1]] <- 0
orig_df$CD62L[orig_df$CD62L<thresh$CD62L[1]] <- 0

int_cols <- c(which(colnames(orig_df)=="CD44"))


cd19_cd44 <- plotter2D(orig_df, "umap", "umap_1", "umap_2", 
                     int_cols, umap_xlims, umap_ylims, splitter=T)

input_df=orig_df
new_df = data.frame(Environment=NA, perc_cd44=NA)
for(k in 1:length(unique(input_df$id))){
  curr <- subset(input_df, id == unique(input_df$id)[k])
  up = length(which(curr[,int_cols]>0))/nrow(curr)*100
  adder <- data.frame(Environment=curr$Environment[1], perc_cd44=up)
  new_df <- rbind(new_df, adder)
}
new_df <- new_df[-1,]

cd19_cd44_box = ggplot(new_df, aes(Environment, perc_cd44, color=Environment)) +
  geom_boxplot(alpha=0.2, outlier.shape = NA) + geom_jitter(width=0.2) +
  scale_color_manual(values=c("mediumorchid3", "red3"))+
  ylab("% of CD44+ cells per mouse") + xlab("")+
  theme_bw() +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12, color='black'),
        #legend.title = element_text(size=15),
        legend.position = 'none',
        legend.title = element_blank(),
        legend.text = element_text(size=12))
 
 #
#
#

#
reader4 <- paste0("inputs/CD4_umap_combo_Blood.txt")
umap_blood_df <- read.table(reader4, F, '\t')
meta <- read.table("mice_metadata.11.19_mouse_id.txt", T, '\t')
reader4 <- paste0("inputs/CD4_Blood_df.txt")
orig_df <- read.table(reader4, F, '\t')
ids <- orig_df$V1

names <- colnames(read.table("name_change.csv", T, ","))

colnames(orig_df) <- c("id", names)

### add metadata

orig_df$id <- factor(orig_df$id, levels = unique(orig_df$id))
orig_df <- orig_df[order(orig_df$id),]

rownames(meta) <- meta$mouse_id
meta <- meta[levels(orig_df$id),]

uniq_ids <- unique(orig_df$id)
orig_df2 <- data.frame(orig_df[1,], Genotype = NA, Environment = NA,
                       Wedge_cage = NA, Gender = NA, Pregnant = NA, Diarrhea = NA, Flow.date=NA)
for (j in 1:length(uniq_ids)){
  curr <- subset(orig_df, id == uniq_ids[j])
  meta_curr <- subset(meta, mouse_id == as.character(uniq_ids[j]))
  
  curr$Genotype <- rep(meta_curr$Genotype, each = nrow(curr))
  curr$Environment <- rep(meta_curr$Environment, each =nrow(curr))
  curr$Wedge_cage <- rep(meta_curr$Wedge_cage, each =nrow(curr))
  curr$Gender <- rep(meta_curr$Gender, each =nrow(curr))
  curr$Pregnant <- rep(meta_curr$Pregnant, each =nrow(curr))
  curr$Diarrhea <- rep(meta_curr$Diarrhea, each =nrow(curr))
  curr$Flow.date <- rep(meta_curr$Flow.date, each =nrow(curr))
  
  orig_df2 <- rbind(orig_df2, curr)
}
orig_df <- orig_df2[-1,]

colnames(orig_df)

# blood all

orig_df$umap_1 <- umap_blood_df$V1
orig_df$umap_2 <- umap_blood_df$V2

umap_xlims <- c(floor(min(orig_df$umap_1)), ceiling(max(orig_df$umap_1)))
umap_ylims <- c(floor(min(orig_df$umap_2)), ceiling(max(orig_df$umap_2)))

reader <- paste0("Distributions/thresholds/CD4_Blood_minor_thresholds.txt")
thresh <- read.table(reader, T, "\t")

orig_df$CD43_1B11[orig_df$CD43_1B11<thresh$CD43_1B11[1]] <- 0
orig_df$PD1[orig_df$PD1<thresh$PD1[1]] <- 0
orig_df$CD25[orig_df$CD25<thresh$CD25[1]] <- 0
orig_df$CD44[orig_df$CD44<thresh$CD44[1]] <- 0
orig_df$CD127[orig_df$CD127<thresh$CD127[1]] <- 0
orig_df$CXCR3[orig_df$CXCR3<thresh$CXCR3[1]] <- 0
orig_df$KLRG1[orig_df$KLRG1<thresh$KLRG1[1]] <- 0
orig_df$CD27[orig_df$CD27<thresh$CD27[1]] <- 0
orig_df$CD69[orig_df$CD69<thresh$CD69[1]] <- 0
orig_df$CD62L[orig_df$CD62L<thresh$CD62L[1]] <- 0

int_cols <- c(which(colnames(orig_df)=="CD62L"))


cd4_cd62L <- plotter2D(orig_df, "umap", "umap_1", "umap_2", 
                       int_cols, umap_xlims, umap_ylims, splitter=T)

input_df=orig_df
new_df = data.frame(Environment=NA, perc_cd44=NA)
for(k in 1:length(unique(input_df$id))){
  curr <- subset(input_df, id == unique(input_df$id)[k])
  up = length(which(curr[,int_cols]>0))/nrow(curr)*100
  adder <- data.frame(Environment=curr$Environment[1], perc_cd44=up)
  new_df <- rbind(new_df, adder)
}
new_df <- new_df[-1,]

cd4_cd62L_box = ggplot(new_df, aes(Environment, perc_cd44, color=Environment)) +
  geom_boxplot(alpha=0.2, outlier.shape = NA) + geom_jitter(width=0.2) +
  scale_color_manual(values=c("mediumorchid3", "red3"))+
  ylab("% of CD62L+ cells per mouse") + xlab("")+
  theme_bw() +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12, color='black'),
        #legend.title = element_text(size=15),
        legend.position = 'none',
        legend.title = element_blank(),
        legend.text = element_text(size=12))

#
#
g1=ggplot()+theme_void()
## boxplots and smaller umaps

p1 <- arrangeGrob(cd19, cd19_cd44_box, nrow=2)
p2 <- arrangeGrob(cd4, cd4_cd62L_box, nrow=2)

p1_void <- arrangeGrob(cd19, g1, nrow=2)
p2_void <- arrangeGrob(cd4, g1, nrow=2)

##



png("plots/Fig1b_c.png", 
    height = 8, width = 14, units = 'in', res=300)
grid.arrange(p1, g1, cd19_cd44[[1]], g1, p2, g1, cd4_cd62L[[1]],
             nrow = 1, widths=c(1,0.2,1,0.2,1,0.2,1))
dev.off()

png("plots/Fig1b_c_void.png", 
    height = 8, width = 14, units = 'in', res=300)
grid.arrange(p1_void, g1, cd19_cd44[[1]], g1, p2_void, g1, cd4_cd62L[[1]],
             nrow = 1, widths=c(1,0.2,1,0.2,1,0.2,1))
dev.off()


#TOO BIG
#pdf("plots/Fig1b_c.pdf", 
#    height = 8, width = 14)
#grid.arrange(p1, g1, cd19_cd44[[1]], g1, p2, g1, cd4_cd62L[[1]],
#             nrow = 1, widths=c(1,0.2,1,0.2,1,0.2,1))
#dev.off()

# png("plots/Fig1b_c_1.png", 
#     height = 5, width = 5, units = 'in', res=500)
# grid.arrange(cd19)
# dev.off()
# 
# 
pdf("plots/Fig1b_c_2.pdf",
    height = 5, width = 5)
grid.arrange(cd19_cd44_box)
dev.off()
# 
# png("plots/Fig1b_c_3.png", 
#     height = 10, width = 5, units = 'in', res=500)
# grid.arrange(cd19_cd44[[1]])
# dev.off()
# 
# 
# png("plots/Fig1b_c_4.png", 
#     height = 5, width = 5, units = 'in', res=500)
# grid.arrange(cd4)
# dev.off()
# 
# png("plots/Fig1b_c_5.png", 
#     height = 5, width = 5, units = 'in', res=500)
# grid.arrange(cd4_cd62L_box)
# dev.off()
# 
pdf("plots/Fig1b_c_5.pdf",
    height = 5, width = 5)
grid.arrange(cd4_cd62L_box)
dev.off()
# 
# png("plots/Fig1b_c_6.png", 
#     height = 10, width = 5, units = 'in', res=500)
# grid.arrange(cd4_cd62L[[1]])
# dev.off()



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
