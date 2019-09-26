### plot UMAPs for all conditions, figures laid out

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
      ramper = colorRampPalette(c("grey99", "purple4"))(2)
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
tsne_blood_df <- read.table("inputs/ftsne_combo_Blood.txt", F, '\t')
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

orig_df$tsne_1 <- tsne_blood_df$V1
orig_df$tsne_2 <- tsne_blood_df$V2

umap_xlims <- c(floor(min(orig_df$umap_1)), ceiling(max(orig_df$umap_1)))
umap_ylims <- c(floor(min(orig_df$umap_2)), ceiling(max(orig_df$umap_2)))
tsne_xlims <- c(floor(min(orig_df$tsne_1)), ceiling(max(orig_df$tsne_1)))
tsne_ylims <- c(floor(min(orig_df$tsne_2)), ceiling(max(orig_df$tsne_1)))


thresh <- read.table("Distributions/thresholds/Blood_major_thresholds.txt", T, "\t")
orig_df$CD3[orig_df$CD3<thresh$CD3[1]] <- 0
orig_df$CD19[orig_df$CD19<thresh$CD19[1]] <- 0
orig_df$CD4[orig_df$CD4<thresh$CD4[1]] <- 0
orig_df$CD8[orig_df$CD8<thresh$CD8[1]] <- 0

int_cols <- c(which(colnames(orig_df)=="CD19"),
              which(colnames(orig_df)=="CD3"),
              which(colnames(orig_df)=="CD4"),
              which(colnames(orig_df)=="CD8"))

cd19_3_4_8_blood <- plotter2D(orig_df, "umap", "umap_1", "umap_2", 
                              int_cols, umap_xlims, umap_ylims, splitter=F)

row1_right <- arrangeGrob(grobs=cd19_3_4_8_blood, nrow=1)

ramper_1 = colorRampPalette(c("purple4"))(1)
orig_df$UMAP <- rep(1, nrow(orig_df))
int_cols <- which(colnames(orig_df)=="UMAP")
umap_all <- plotter2D(orig_df, "umap", "umap_1", "umap_2", 
                      int_cols, umap_xlims, umap_ylims, splitter=F, ramper=ramper_1)
row1_umap <- arrangeGrob(grobs=umap_all)

orig_df$tSNE <- rep(1, nrow(orig_df))
int_cols <- which(colnames(orig_df)=="tSNE")
tsne_all <- plotter2D(orig_df, "tsne", "tsne_1", "tsne_2", 
                      int_cols, tsne_xlims, tsne_ylims, splitter=F, ramper=ramper_1)
row1_tsne <- arrangeGrob(grobs=tsne_all)

row1 <- arrangeGrob(row1_tsne, row1_umap, row1_right, nrow=1, widths = c(1,1,4))

#
#
#
#
#

#
#
#
### Blood CD19, CD4, CD8

marks <- c("CD19", "CD4", "CD8")
row_adder = list()
for(i in 1:length(marks)){
  
  reader <- paste0("inputs/", marks[i], "_umap_combo_Blood.txt")
  umap_blood_df <- read.table(reader, F, '\t')
  meta <- read.table("mice_metadata.11.19_mouse_id.txt", T, '\t')
  reader <- paste0("inputs/", marks[i], "_Blood_df.txt")
  orig_df <- read.table(reader, F, '\t')
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
  
  reader <- paste0("Distributions/thresholds/", marks[i], "_Blood_minor_thresholds.txt")
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
  
  int_cols <- c(which(colnames(orig_df)=="CD43_1B11"),
                which(colnames(orig_df)=="PD1"),
                which(colnames(orig_df)=="CD25"),
                which(colnames(orig_df)=="CD44"),
                which(colnames(orig_df)=="CD127"),
                which(colnames(orig_df)=="CXCR3"),
                which(colnames(orig_df)=="KLRG1"),
                which(colnames(orig_df)=="CD27"),
                which(colnames(orig_df)=="CD69"),
                which(colnames(orig_df)=="CD62L"))
  
  namer <- paste0("Distributions/", marks[i], "_Blood_distr.pdf")
  pdf(namer, height=8,width=15)
  par(mfrow = c(2, 5))
  for(j in int_cols){
    
    ylims=c(0,0.01)
    plot(density(orig_df[,j]), 
         main=paste0(colnames(orig_df)[j]),ylim=ylims)
  }
  dev.off()
  
  new_row <- plotter2D(orig_df, "umap", "umap_1", "umap_2", 
                                int_cols, umap_xlims, umap_ylims, splitter=T)
  
  row_adder[[i]] <- arrangeGrob(grobs=new_row, nrow=1)
}

row2 <- row_adder[[1]]
row3 <- row_adder[[2]]
row4 <- row_adder[[3]]

png("plots/Blood_Lymph_all.png", 
    height = 11, width = 11, units = 'in', res=400)
grid.arrange(row1, row2, row3, row4,
             nrow = 4, heights=c(1,2,2,2))



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
# MLN Lymph
# Row 1. TSNE vs UMAP & CD3, CD19, CD4, CD8 Combined

### MLN all

umap_blood_df <- read.table("inputs/umap_combo_MLN.txt", F, '\t')
tsne_blood_df <- read.table("inputs/ftsne_combo_MLN.txt", F, '\t')
meta <- read.table("mice_metadata.11.19_mouse_id.txt", T, '\t')
orig_df <- read.table("inputs/MLN_df.txt", F, '\t')
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

orig_df$tsne_1 <- tsne_blood_df$V1
orig_df$tsne_2 <- tsne_blood_df$V2

umap_xlims <- c(floor(min(orig_df$umap_1)), ceiling(max(orig_df$umap_1)))
umap_ylims <- c(floor(min(orig_df$umap_2)), ceiling(max(orig_df$umap_2)))
tsne_xlims <- c(floor(min(orig_df$tsne_1)), ceiling(max(orig_df$tsne_1)))
tsne_ylims <- c(floor(min(orig_df$tsne_2)), ceiling(max(orig_df$tsne_1)))


thresh <- read.table("Distributions/thresholds/MLN_major_thresholds.txt", T, "\t")
orig_df$CD3[orig_df$CD3<thresh$CD3[1]] <- 0
orig_df$CD19[orig_df$CD19<thresh$CD19[1]] <- 0
orig_df$CD4[orig_df$CD4<thresh$CD4[1]] <- 0
orig_df$CD8[orig_df$CD8<thresh$CD8[1]] <- 0

int_cols <- c(which(colnames(orig_df)=="CD19"),
              which(colnames(orig_df)=="CD3"),
              which(colnames(orig_df)=="CD4"),
              which(colnames(orig_df)=="CD8"))

cd19_3_4_8_blood <- plotter2D(orig_df, "umap", "umap_1", "umap_2", 
                              int_cols, umap_xlims, umap_ylims, splitter=F)

row1_right <- arrangeGrob(grobs=cd19_3_4_8_blood, nrow=1)

ramper_1 = colorRampPalette(c("purple4"))(1)
orig_df$UMAP <- rep(1, nrow(orig_df))
int_cols <- which(colnames(orig_df)=="UMAP")
umap_all <- plotter2D(orig_df, "umap", "umap_1", "umap_2", 
                      int_cols, umap_xlims, umap_ylims, splitter=F, ramper=ramper_1)
row1_umap <- arrangeGrob(grobs=umap_all)

orig_df$tSNE <- rep(1, nrow(orig_df))
int_cols <- which(colnames(orig_df)=="tSNE")
tsne_all <- plotter2D(orig_df, "tsne", "tsne_1", "tsne_2", 
                      int_cols, tsne_xlims, tsne_ylims, splitter=F, ramper=ramper_1)
row1_tsne <- arrangeGrob(grobs=tsne_all)

row1 <- arrangeGrob(row1_tsne, row1_umap, row1_right, nrow=1, widths = c(1,1,4))

#
#
#
#
#

#
#
#
### Blood CD19, CD4, CD8

marks <- c("CD19", "CD4", "CD8")
row_adder = list()
for(i in 1:length(marks)){
  
  reader <- paste0("inputs/", marks[i], "_umap_combo_MLN.txt")
  umap_blood_df <- read.table(reader, F, '\t')
  meta <- read.table("mice_metadata.11.19_mouse_id.txt", T, '\t')
  reader <- paste0("inputs/", marks[i], "_MLN_df.txt")
  orig_df <- read.table(reader, F, '\t')
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
  
  reader <- paste0("Distributions/thresholds/", marks[i], "_MLN_minor_thresholds.txt")
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
  
  int_cols <- c(which(colnames(orig_df)=="CD43_1B11"),
                which(colnames(orig_df)=="PD1"),
                which(colnames(orig_df)=="CD25"),
                which(colnames(orig_df)=="CD44"),
                which(colnames(orig_df)=="CD127"),
                which(colnames(orig_df)=="CXCR3"),
                which(colnames(orig_df)=="KLRG1"),
                which(colnames(orig_df)=="CD27"),
                which(colnames(orig_df)=="CD69"),
                which(colnames(orig_df)=="CD62L"))
  
  namer <- paste0("Distributions/", marks[i], "_MLN_distr.pdf")
  pdf(namer, height=8,width=15)
  par(mfrow = c(2, 5))
  for(j in int_cols){
    
    ylims=c(0,0.01)
    plot(density(orig_df[,j]), 
         main=paste0(colnames(orig_df)[j]),ylim=ylims)
  }
  dev.off()
  
  new_row <- plotter2D(orig_df, "umap", "umap_1", "umap_2", 
                       int_cols, umap_xlims, umap_ylims, splitter=T)
  
  row_adder[[i]] <- arrangeGrob(grobs=new_row, nrow=1)
}

row2 <- row_adder[[1]]
row3 <- row_adder[[2]]
row4 <- row_adder[[3]]

png("plots/MLN_Lymph_all.png", 
    height = 11, width = 11, units = 'in', res=400)
grid.arrange(row1, row2, row3, row4,
             nrow = 4, heights=c(1,2,2,2))



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
#
# MLN Myeloid
# Row 1. TSNE vs UMAP & CD3, CD19, CD4, CD8 Combined

### MLN myeloid

# umap_blood_df <- read.table("inputs/umap_combo_MLN_myeloid.txt", F, '\t')
# tsne_blood_df <- read.table("inputs/ftsne_combo_MLN_myeloid.txt", F, '\t')
# meta <- read.table("mice_metadata.11.19_mouse_id.txt", T, '\t')
# orig_df <- read.table("inputs/MLN_myeloid_df.txt", F, '\t')
# ids <- orig_df$V1
# 
# names <- colnames(read.table("name_change_myeloid.csv", T, ","))
# 
# colnames(orig_df) <- c("id", names)
# 
# ### add metadata
# 
# orig_df$id <- factor(orig_df$id, levels = unique(orig_df$id))
# orig_df <- orig_df[order(orig_df$id),]
# 
# rownames(meta) <- meta$mouse_id
# meta <- meta[levels(orig_df$id),]
# 
# uniq_ids <- unique(orig_df$id)
# orig_df2 <- data.frame(orig_df[1,], Genotype = NA, Environment = NA,
#                        Wedge_cage = NA, Gender = NA, Pregnant = NA, Diarrhea = NA, Flow.date=NA)
# for (j in 1:length(uniq_ids)){
#   curr <- subset(orig_df, id == uniq_ids[j])
#   meta_curr <- subset(meta, mouse_id == as.character(uniq_ids[j]))
#   
#   curr$Genotype <- rep(meta_curr$Genotype, each = nrow(curr))
#   curr$Environment <- rep(meta_curr$Environment, each =nrow(curr))
#   curr$Wedge_cage <- rep(meta_curr$Wedge_cage, each =nrow(curr))
#   curr$Gender <- rep(meta_curr$Gender, each =nrow(curr))
#   curr$Pregnant <- rep(meta_curr$Pregnant, each =nrow(curr))
#   curr$Diarrhea <- rep(meta_curr$Diarrhea, each =nrow(curr))
#   curr$Flow.date <- rep(meta_curr$Flow.date, each =nrow(curr))
#   
#   orig_df2 <- rbind(orig_df2, curr)
# }
# orig_df <- orig_df2[-1,]
# 
# colnames(orig_df)
# 
# # blood all
# 
# orig_df$umap_1 <- umap_blood_df$V1
# orig_df$umap_2 <- umap_blood_df$V2
# 
# orig_df$tsne_1 <- tsne_blood_df$V1
# orig_df$tsne_2 <- tsne_blood_df$V2
# 
# umap_xlims <- c(floor(min(orig_df$umap_1)), ceiling(max(orig_df$umap_1)))
# umap_ylims <- c(floor(min(orig_df$umap_2)), ceiling(max(orig_df$umap_2)))
# tsne_xlims <- c(floor(min(orig_df$tsne_1)), ceiling(max(orig_df$tsne_1)))
# tsne_ylims <- c(floor(min(orig_df$tsne_2)), ceiling(max(orig_df$tsne_1)))
# 
# 
# #thresh <- read.table("Distributions/thresholds/MLN_myeloid_major_thresholds.txt", T, "\t")
# #orig_df$CD3[orig_df$CD3<thresh$CD3[1]] <- 0
# #orig_df$CD19[orig_df$CD19<thresh$CD19[1]] <- 0
# #orig_df$CD4[orig_df$CD4<thresh$CD4[1]] <- 0
# #orig_df$CD8[orig_df$CD8<thresh$CD8[1]] <- 0
# 
# # int_cols <- c(which(colnames(orig_df)=="CD19"),
# #               which(colnames(orig_df)=="CD3"),
# #               which(colnames(orig_df)=="CD4"),
# #               which(colnames(orig_df)=="CD8"))
# 
# int_cols <- which(colnames(orig_df) %in% c("CD40","Ly6G","Ly6C",
#                                            "PDL1","CD64","CD11c",
#                                            "CD103","CD69","Siglec_F",
#                                            "CD86",	"CD11b", "F4_80", 
#                                            "MHCII","PDL2")
# )
# 
# namer <- paste0("Distributions/MLN_myeloid_distr.pdf")
# pdf(namer, height=8,width=20)
# par(mfrow = c(2, 7))
# for(j in int_cols){
#   
#   ylims=c(0,0.01)
#   plot(density(orig_df[,j]), 
#        main=paste0(colnames(orig_df)[j]),ylim=ylims)
# }
# dev.off()
# 
# all_myeloid <- plotter2D(orig_df, "umap", "umap_1", "umap_2", 
#                               int_cols, umap_xlims, umap_ylims, splitter=F, ramper=ramper_more)
# 
# row1_right <- arrangeGrob(grobs=all_myeloid, nrow=2)
# 
# ramper_1 = colorRampPalette(c("purple4"))(1)
# orig_df$UMAP <- rep(1, nrow(orig_df))
# int_cols <- which(colnames(orig_df)=="UMAP")
# umap_all <- plotter2D(orig_df, "umap", "umap_1", "umap_2", 
#                       int_cols, umap_xlims, umap_ylims, splitter=F, ramper=ramper_1)
# row1_umap <- arrangeGrob(grobs=umap_all)
# 
# orig_df$tSNE <- rep(1, nrow(orig_df))
# int_cols <- which(colnames(orig_df)=="tSNE")
# tsne_all <- plotter2D(orig_df, "tsne", "tsne_1", "tsne_2", 
#                       int_cols, tsne_xlims, tsne_ylims, splitter=F, ramper=ramper_1)
# row1_tsne <- arrangeGrob(grobs=tsne_all)
# 
# #row1 <- arrangeGrob(row1_tsne, row1_umap, row1_right, nrow=1, widths = c(1,1,4))
# 
# row_left <- arrangeGrob(row1_tsne, row1_umap, nrow=2)
# 
# png("plots/MLN_myeloid_all.png", 
#     height = 5.5, width = 18, units = 'in', res=400)
# grid.arrange(row_left, row1_right,
#              nrow = 1, widths=c(1,7))
# dev.off()


