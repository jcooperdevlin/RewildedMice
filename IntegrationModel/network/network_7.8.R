#### network creator

setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/IntegrationModel/network")

library(igraph)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(gridExtra)
library(ggplotify)
library(scales)
library(ggsignif)

genes = readRDS("../ml_inputs/gene_modules.RData")
cytokines = readRDS("../ml_inputs/cytokines.RData")
cells = readRDS("../ml_inputs/facs.RData")
#otus = readRDS("../ml_inputs/otus_relab.RData")
otus = readRDS("../ml_inputs/otus.RData")

p_mats <- list.files("sim_mats", full.names = T, pattern = c("p.adj"))
c_mats <- setdiff(list.files("sim_mats", pattern = c("txt"), full.names = T),p_mats)

source("sim_mat_func.R")

mats = sim_combine(c_mats, p_mats, 0.01)



#

#
#
#
#
net_mat <- mats[[1]]

##filter
net_mat[net_mat<0.6 & net_mat > -0.6]=0
good_cols <- which(colSums(abs(net_mat))>0)
net_mat=net_mat[good_cols,good_cols]

feats_df <- mats[[2]][good_cols,]

#store 4 later 
net_mat_og <- net_mat

#otu_cols <- grep("D_0__Bacteria", colnames(net_mat))
#colnames(net_mat)[otu_cols] <- gsub(".*__", "", colnames(net_mat)[otu_cols])

#make facs names shorter
#mln_cols <- grep("CD44", colnames(net_mat))
#colnames(net_mat)[mln_cols] <- gsub("T_cells_", "T_cells\n", colnames(net_mat)[mln_cols])

#mln_cols <- grep("MHCII", colnames(net_mat))
#colnames(net_mat)[mln_cols] <- gsub("CD11c_", "CD11c\n", colnames(net_mat)[mln_cols])

#
#

#



##

pp = prcomp(net_mat, scale = F)
p_plot <- data.frame(PC1=pp$x[,1], PC2=pp$x[,2], feat=colnames(net_mat))
ggplot(p_plot, aes(PC1, PC2, label=feat)) + geom_label()

connecs = c()
for(i in 1:ncol(net_mat)){
  huh = length(which(net_mat[,i]!=0))
  connecs <- c(connecs, huh)
}

p_plot$edges <- connecs

ggplot(p_plot, aes(PC1, PC2, label=edges)) + geom_label()

ggplot(p_plot, aes(PC1, PC2, size=edges)) + geom_point() +
  scale_radius(breaks = seq(1,35,4))

p_plot$new_x <- p_plot$PC1
p_plot$new_y <- p_plot$PC2
for(i in 1:nrow(p_plot)){
  use_plot <- p_plot[-i,]
  
  curr_x <- p_plot$PC1[i]
  curr_y <- p_plot$PC2[i]
  
  up_t1 <- curr_x+0.5
  down_t1 <- curr_x-0.5
  up_t2 <- curr_y+0.5
  down_t2 <- curr_y-0.5
  
  close_point <- subset(use_plot, PC1 < up_t1 & PC1 > down_t1 &
                      PC2 < up_t2 & PC2 > down_t2)
  
  if(nrow(close_point)>0){
    if(curr_x > 0) {
      dirr_x = 1
    } else {dirr_x = -1}
    if(curr_y > 0){
      dirr_y = 1
    } else {dirr_y = -1}
    
    adder_x <- dirr_x*runif(1, -1, 1)
    adder_y <- dirr_y*runif(1, -1, 1)
    p_plot$new_x[i] <- p_plot$new_x[i]+adder_x
    p_plot$new_y[i] <- p_plot$new_y[i]+adder_y
    
  }
}

ggplot(p_plot, aes(new_x, new_y, size=edges)) + geom_point() +
  scale_radius(breaks = seq(1,35,4))


for(i in 1:nrow(p_plot)){
  use_plot <- p_plot[-i,]
  
  curr_x <- p_plot$new_x[i]
  curr_y <- p_plot$new_y[i]
  
  up_t1 <- curr_x+0.5
  down_t1 <- curr_x-0.5
  up_t2 <- curr_y+0.5
  down_t2 <- curr_y-0.5
  
  close_point <- subset(use_plot, PC1 < up_t1 & PC1 > down_t1 &
                          PC2 < up_t2 & PC2 > down_t2)
  
  if(nrow(close_point)>0){
    if(curr_x > 0) {
      dirr_x = 1
    } else {dirr_x = -1}
    if(curr_y > 0){
      dirr_y = 1
    } else {dirr_y = -1}
    
    adder_x <- dirr_x*runif(1, 0, 1)
    adder_y <- dirr_y*runif(1, 0, 1)
    p_plot$new_x[i] <- p_plot$new_x[i]+adder_x
    p_plot$new_y[i] <- p_plot$new_y[i]+adder_y
    
  }
}

ggplot(p_plot, aes(new_x, new_y, size=edges)) + geom_point() +
  scale_radius(breaks = seq(1,35,4))

p_plot$edge_adj <- (round(p_plot$edges/10,0)+1)*3

p_plot$new_name <- as.character(p_plot$feat)
p_plot$new_name[p_plot$edges<3] <- ""
### think about ways to repel data points but keep PCA structure

#### READ IN GOOD COORDINATES!!!!
good_net_coords <- read.table("network7.8_omics_coords.txt", header=T, sep='\t')
p_plot <- good_net_coords

#
#
#
#
#

network=graph_from_adjacency_matrix(as.matrix(net_mat), weighted=T, mode="undirected", diag=F)
icoord = layout_with_fr(network)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

E(network)$color[E(network)$weight > 0] <- add.alpha('red', alpha=0.25)
E(network)$color[E(network)$weight < 0] <- add.alpha('blue', alpha=0.25)


###calculate edge thickness
E(network)

edge_df <- data.frame(as_edgelist(network), weight = E(network)$weight)
edge_df$thickness <- rescale(abs(edge_df$weight), to=c(0,8))

#
#
#
#


network <- set_vertex_attr(network, "label", value=as.character(p_plot$new_name))

coul <- c(cytokines="darkseagreen2", genes="lightcoral", cells="dodgerblue2", otus="mediumorchid3")
my_color=coul[feats_df$dataType]
# plot

pdf("full_network_by_pca_7.8.pdf", height = 50, width = 50)
par(bg="white", mar=c(0,0,0,0))
set.seed(4)
plot(network, #layout=icoord,
     layout=as.matrix(p_plot[,5:6]),
     #vertex.size=5,
     vertex.size=as.matrix(p_plot$edge_adj),
     vertex.shape='circle',
     vertex.color=my_color, 
     vertex.label.cex=1.2,
     vertex.label.family="Helvetica",
     vertex.label.color="black",
     vertex.frame.color="transparent"
)
dev.off()


curr_coord = icoord
new_coord = data.frame(curr_coord, feats_df)
colnames(new_coord) <- c("X", "Y", 'Feature', "DataType")

library(ggplot2)
ggplot(new_coord, aes(X,Y,color=DataType)) +
  geom_point()

gs = as.character(unique(new_coord$DataType))
gs=c("genes", "cells", "cytokines", "otus")
better_coords <- data.frame(X=NA,Y=NA,Feature=NA,DataType=NA,X2=NA,Y2=NA)
low_bounds_y=c(250,800,250,0)
high_bounds_y=c(700,950,700,200)

low_bounds_x=c(0,300,700,300)
high_bounds_x=c(300,700,1000,700)

for(i in 1:length(gs)){
  curr <- subset(new_coord, DataType == gs[i])

  x2 = sample(low_bounds_x[i]:high_bounds_x[i], nrow(curr))
  y2 = sample(low_bounds_y[i]:high_bounds_y[i], nrow(curr), replace = T)

  curr$X2=x2
  curr$Y2=y2

  better_coords = rbind(better_coords, curr)
}
better_coords=better_coords[-1,]

ggplot(better_coords, aes(X2,Y2,color=DataType)) +
  geom_point()
better_coords=better_coords[rownames(new_coord),]

pdf("full_network_omic_7.8.pdf", height = 50, width = 50)
par(bg="white", mar=c(0,0,0,0))
set.seed(4)
plot(network, layout=as.matrix(better_coords[,5:6]),
     #vertex.size=10,
     vertex.size=as.matrix(p_plot$edge_adj),
     vertex.color=my_color,
     vertex.label.cex=1.2,
     vertex.label.family="Helvetica",
     vertex.label.color="black",
     vertex.frame.color="transparent"
)
dev.off()

tkplot(network, #layout=as.matrix(better_coords[,5:6]),
     layout=as.matrix(p_plot[,9:10]),
     #vertex.size=5,
     vertex.size=as.matrix(p_plot$edge_adj),
     vertex.shape='circle',
     vertex.color=my_color, 
     vertex.label.cex=1.2,
     vertex.label.family="Helvetica",
     vertex.label.color="black"#,
     #vertex.frame.color="transparent"
)

#saveRDS(network, "FINAL_NETWORK_7.22.Rdata")
#readRDS("FINAL_NETWORK_7.22.Rdata")

tk_coords = tk_coords(3)

pdf("full_network_omic_7.8.pdf", height = 50, width = 50)
par(bg="white", mar=c(0,0,0,0))
set.seed(4)
plot(network, #layout=as.matrix(better_coords[,5:6]),
     #layout=as.matrix(p_plot[,5:6]),
     #vertex.size=5,
     layout=tk_coords,
     edge.width = as.matrix(edge_df$thickness),
     vertex.size=as.matrix(p_plot$edge_adj),
     vertex.shape='circle',
     vertex.color=my_color, 
     vertex.label.cex=1.2,
     vertex.label.family="Helvetica",
     vertex.label.color="black",
     vertex.frame.color="transparent"
)

#legend(x=-1.5, y=-0.6, 
#       legend=paste( levels(as.factor(feats_df$DataType[good_cols]))), 
#       col = coul , bty = "n", pch=20 , pt.cex = 14, cex = 9, text.col="black" , horiz = F)
dev.off()



#
#
#
#### SAVE THE COORDINATES

##as of 7.8 here is the network!!!!
p_plot <- cbind(p_plot, tk_coords)
write.table(p_plot, "network7.8_omics_coords.txt", quote=F, sep='\t', row.names=F)
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
### stats about this network

nums <- unlist(c(net_mat))
rows <- rep(rownames(net_mat), ncol(net_mat))
cols <- rep(colnames(net_mat), each=nrow(net_mat))
type_row <- rep(feats_df$dataType, ncol(net_mat))
type_col <- rep(feats_df$dataType, each=nrow(net_mat))

int_df <- data.frame(rows, cols, type_row, type_col, nums)
int_real <- subset(int_df, nums != 0)
table(int_real$type_row, int_real$type_col)

popular_df <- data.frame(table(int_real$rows))

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
#
#
#
#
#

#
###
# GSA
genesets=list()
num=1
for(i in 1:nrow(net_mat)){
  curr = net_mat[i,]
  gs=c(rownames(net_mat)[i],rownames(net_mat)[which(curr != 0)])
  
  if(length(gs)>3){
    print(length(gs))
    genesets[[num]]<-gs
    num=num+1
  }
}



gnames=paste0("module_", 1:length(genesets))

custom_genesets = list(genesets = genesets, geneset.names = gnames, 
                       geneset.descriptions = gnames)
class(custom_genesets) = "GSA.genesets"

library(GSA)
set.seed(100)

colnames(otus)=make.names(colnames(otus))
x=data.frame(rbind(t(genes), t(cytokines), t(cells), t(otus)))
colnames(x) <- gsub("X", "", colnames(x))

meta <- readRDS("../ml_inputs/mouse_metadata.RData")
y=as.character(meta$Environment)
y[y=="lab"]<-1
y[y=="wild"]<-2
y=as.numeric(y)

GSA.obj<-GSA(x,y, genenames = rownames(x), genesets=custom_genesets[[1]], 
             resp.type="Two class unpaired", nperms=100, minsize = 10)
GSA.listsets(GSA.obj, geneset.names=gnames,FDRcut=0.1)
GSA.plot(GSA.obj)
good_gs=rbind(GSA.listsets(GSA.obj, geneset.names=gnames,FDRcut=0.1)$positive, 
          GSA.listsets(GSA.obj, geneset.names=gnames,FDRcut=0.1)$negative)

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
###### given list of sig gene sets make a subnetwork and some boxplot comparisons?

genesets <- good_gs[,1]
GSA_networker <- function(genesets, x_data, y_data, meta_entry, box_cat, box_cols){
  glist=list()
  for(i in 1:length(genesets)){
  good_genes <- custom_genesets[[1]][[as.numeric(genesets[i])]]
  
  new_ids <- rownames(net_mat_og)[which(rownames(net_mat_og)%in%good_genes)]
  
  new_net <- net_mat_og[new_ids,new_ids]
  up_feat <- mats[[2]][new_ids,]
  
  network=graph_from_adjacency_matrix(as.matrix(new_net), weighted=T, mode="undirected", diag=F)
  
  icoord = layout_with_fr(network)
  
  add.alpha <- function(col, alpha=1){
    if(missing(col))
      stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2, 
          function(x) 
            rgb(x[1], x[2], x[3], alpha=alpha))  
  }
  
  E(network)$color[E(network)$weight > 0] <- add.alpha('red', alpha=0.2)
  E(network)$color[E(network)$weight < 0] <- add.alpha('blue', alpha=0.2)
  
  coul <- c(cytokines="darkseagreen2", genes="lightcoral", cells="dodgerblue2", otus="mediumorchid3")
  my_color=coul[up_feat$dataType]
  # plot
  
  curr_coord = icoord
  new_coord = data.frame(curr_coord, up_feat)
  colnames(new_coord) <- c("X", "Y", 'Feature', "DataType")
  
  #library(ggplot2)
  #ggplot(new_coord, aes(X,Y,color=DataType)) + 
  #  geom_point()
  
  gs = as.character(unique(new_coord$DataType))
  gs=c("genes", "cells", "cytokines", "otus")
  better_coords <- data.frame(X=NA,Y=NA,Feature=NA,DataType=NA,X2=NA,Y2=NA)
  low_bounds_y=c(250,800,250,0)
  high_bounds_y=c(700,950,700,200)
  
  low_bounds_x=c(0,300,700,300)
  high_bounds_x=c(300,700,1000,700)
  
  for(j in 1:length(gs)){
    curr <- subset(new_coord, DataType == gs[j])
    
    x2 = sample(low_bounds_x[j]:high_bounds_x[j], nrow(curr))
    y2 = sample(low_bounds_y[j]:high_bounds_y[j], nrow(curr), replace = T)
    
    curr$X2=x2
    curr$Y2=y2
    
    better_coords = rbind(better_coords, curr)
  }
  better_coords=better_coords[-1,]
  
  #ggplot(better_coords, aes(X2,Y2,color=DataType)) + 
  #  geom_point()
  better_coords=better_coords[rownames(new_coord),]
  
  network2 = set_edge_attr(network, "label", value="")
  netdraw <- as.grob(function()
    plot(network2, layout=as.matrix(better_coords[,5:6]),
         vertex.size=10,
         vertex.color=my_color, 
         vertex.label.cex=1.1,
         vertex.label.family="Helvetica",
         vertex.label.color="black",
         vertex.frame.color="transparent"
    )
  )
  
  ## and a heat
  
  new_x <- x_data[good_genes,]
  va <- columnAnnotation(df = data.frame(Environment=meta_entry$Environment,
                                         Genotype=meta_entry$Genotype),
                         col = list(Genotype = c(AtgW = "darkorange1",
                                                 AtgE = "dodgerblue2",
                                                 AtgH = "navyblue",
                                                 B6 = "darkorange1",
                                                 NOD2 = "mediumspringgreen"),
                                    Environment=c(wild = "red3", lab = "mediumorchid3")))
  hh=grid.grabExpr(draw(
    Heatmap(t(scale(t(new_x))), top_annotation = va,
             col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
             cluster_columns = T, cluster_rows = T,
             show_row_names = T, show_column_names = F,
             row_names_gp = gpar(fontsize = 10),
             heatmap_legend_param=list(title = ""))
  ))
  
  #
  #
  #
  ##and some boxes
  
  box_df <- data.frame(t(new_x), cat=y_data)
  colnames(box_df) <- c(good_genes, "cat")
  box_df <- melt(box_df)
  sublist=list()
  for(k in 1:length(unique(box_df$variable))){
    plotplot <- subset(box_df, variable == unique(box_df$variable)[k])
    g=ggplot(plotplot, aes(cat, log2(value+1), color=cat, group=cat)) +
      geom_boxplot(alpha=0.5, outlier.shape=NA) +
      geom_jitter(width = 0.2) +
      ggtitle(plotplot$variable[1]) +
      scale_color_manual(values = box_cols) +
      scale_fill_manual(values = box_cols) +
      theme_bw() +
      ylab("Log2 Value") + xlab(box_cat) +
      theme(axis.title = element_text(size=15),
            axis.text = element_text(size=12, color="black"),
            legend.position = 'none',
            legend.title = element_text(size=15),
            legend.text = element_text(size=12, color="black")
      )
    sublist[[k]]<-arrangeGrob(g)
  }
  box_exp <- arrangeGrob(grobs=sublist, nrow=4)
  #box_exp = g+facet_wrap(~variable, ncol = 4)
  
  top_plot <- arrangeGrob(netdraw, hh, nrow=1)
  final_plot <- arrangeGrob(top_plot, box_exp, nrow=2, heights = c(2,3))
  glist[[i]] <- final_plot
  }
  return(glist)
}

hey = GSA_networker(good_gs[,1], x_data=x, y_data = meta$Environment, meta_entry=meta,
                    box_cat = "Environment", 
                    box_cols = c("lab"="mediumorchid3", "wild"="red3"))

pdf("lab_v_wild_nets.pdf", height = 20, width = 20)
marrangeGrob(grobs=hey,ncol=1,nrow=1)
dev.off()



###grab CCL2 candida network and make it nice

custom_genesets[[1]][102]

good_genes <- c("IL.1a_BacteroidesV", custom_genesets[[1]][[102]])

new_ids <- rownames(net_mat_og)[which(rownames(net_mat_og)%in%good_genes)]

new_net <- net_mat_og[new_ids,new_ids]
up_feat <- mats[[2]][new_ids,]

## edge thickness
connecs = c()
for(i in 1:ncol(new_net)){
  huh = length(which(new_net[,i]!=0))
  connecs <- c(connecs, huh)
}

p_plot <- data.frame(names=rownames(new_net))
p_plot$edges = connecs
p_plot$edge_adj <- (round(p_plot$edges/2,0)+1)*3

network=graph_from_adjacency_matrix(as.matrix(new_net), weighted=T, mode="undirected", diag=F)

icoord = layout_with_fr(network)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

E(network)$color[E(network)$weight > 0] <- add.alpha('red', alpha=0.2)
E(network)$color[E(network)$weight < 0] <- add.alpha('blue', alpha=0.2)

coul <- c(cytokines="darkseagreen2", genes="lightcoral", cells="dodgerblue2", otus="mediumorchid3")
my_color=coul[up_feat$dataType]
# plot

curr_coord = icoord
new_coord = data.frame(curr_coord, up_feat)
colnames(new_coord) <- c("X", "Y", 'Feature', "DataType")

#library(ggplot2)
#ggplot(new_coord, aes(X,Y,color=DataType)) + 
#  geom_point()

gs = as.character(unique(new_coord$DataType))
gs=c("genes", "cells", "cytokines", "otus")
better_coords <- data.frame(X=NA,Y=NA,Feature=NA,DataType=NA,X2=NA,Y2=NA)
low_bounds_y=c(250,800,250,0)
high_bounds_y=c(700,950,700,200)

low_bounds_x=c(0,300,700,300)
high_bounds_x=c(300,700,1000,700)

for(j in 1:length(gs)){
  curr <- subset(new_coord, DataType == gs[j])
  
  x2 = sample(low_bounds_x[j]:high_bounds_x[j], nrow(curr))
  y2 = sample(low_bounds_y[j]:high_bounds_y[j], nrow(curr), replace = T)
  
  curr$X2=x2
  curr$Y2=y2
  
  better_coords = rbind(better_coords, curr)
}
better_coords=better_coords[-1,]

#ggplot(better_coords, aes(X2,Y2,color=DataType)) + 
#  geom_point()
better_coords=better_coords[rownames(new_coord),]

network2 = set_edge_attr(network, "label", value="")

tkplot(network2, layout=as.matrix(better_coords[,5:6]),
       vertex.size=p_plot$edge_adj,
       vertex.color=my_color, 
       vertex.label.cex=1.1,
       vertex.label.family="Helvetica",
       vertex.label.color="black"
)

tk_coords <- tk_coords(5)

netdraw <- as.grob(function()
  plot(network2, layout=as.matrix(tk_coords),
       vertex.size=p_plot$edge_adj,
       vertex.color=my_color, 
       vertex.label.cex=1.1,
       vertex.label.family="Helvetica",
       vertex.label.color="black",
       vertex.frame.color="transparent"
  )
)

## and a heat

new_x <- x[good_genes,]
va <- columnAnnotation(df = data.frame(Environment=meta$Environment,
                                       Genotype=meta$Genotype),
                       col = list(Genotype = c(AtgW = "darkorange1",
                                               AtgE = "dodgerblue2",
                                               AtgH = "navyblue",
                                               B6 = "darkorange1",
                                               NOD2 = "mediumspringgreen"),
                                  Environment=c(wild = "red3", lab = "mediumorchid3")))
hh=grid.grabExpr(draw(
  Heatmap(t(scale(t(new_x))), top_annotation = va,
          col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
          cluster_columns = T, cluster_rows = T,
          show_row_names = T, show_column_names = F,
          row_names_gp = gpar(fontsize = 10),
          heatmap_legend_param=list(title = ""))
))

pdf("fig4/cclr_candida_network.pdf", height = 7, width = 20)
grid.arrange(netdraw, hh, nrow=1)
dev.off()



#

#
#
### 11.2 network for genotype specific

custom_genesets[[1]][29]

good_genes <- custom_genesets[[1]][[29]]

new_ids <- rownames(net_mat_og)[which(rownames(net_mat_og)%in%good_genes)]

new_net <- net_mat_og[new_ids,new_ids]
up_feat <- mats[[2]][new_ids,]

network=graph_from_adjacency_matrix(as.matrix(new_net), weighted=T, mode="undirected", diag=F)

icoord = layout_with_fr(network)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

E(network)$color[E(network)$weight > 0] <- add.alpha('red', alpha=0.2)
E(network)$color[E(network)$weight < 0] <- add.alpha('blue', alpha=0.2)

coul <- c(cytokines="darkseagreen2", genes="lightcoral", cells="dodgerblue2", otus="mediumorchid3")
my_color=coul[up_feat$dataType]
# plot

curr_coord = icoord
new_coord = data.frame(curr_coord, up_feat)
colnames(new_coord) <- c("X", "Y", 'Feature', "DataType")

#library(ggplot2)
#ggplot(new_coord, aes(X,Y,color=DataType)) + 
#  geom_point()

gs = as.character(unique(new_coord$DataType))
gs=c("genes", "cells", "cytokines", "otus")
better_coords <- data.frame(X=NA,Y=NA,Feature=NA,DataType=NA,X2=NA,Y2=NA)
low_bounds_y=c(250,800,250,0)
high_bounds_y=c(700,950,700,200)

low_bounds_x=c(0,300,700,300)
high_bounds_x=c(300,700,1000,700)

for(j in 1:length(gs)){
  curr <- subset(new_coord, DataType == gs[j])
  
  x2 = sample(low_bounds_x[j]:high_bounds_x[j], nrow(curr))
  y2 = sample(low_bounds_y[j]:high_bounds_y[j], nrow(curr), replace = T)
  
  curr$X2=x2
  curr$Y2=y2
  
  better_coords = rbind(better_coords, curr)
}
better_coords=better_coords[-1,]

#ggplot(better_coords, aes(X2,Y2,color=DataType)) + 
#  geom_point()
better_coords=better_coords[rownames(new_coord),]

network2 = set_edge_attr(network, "label", value="")

tkplot(network2, layout=as.matrix(better_coords[,5:6]),
       vertex.size=10,
       vertex.color=my_color, 
       vertex.label.cex=1.1,
       vertex.label.family="Helvetica",
       vertex.label.color="black"
)

tk_coords <- tk_coords(5)

netdraw <- as.grob(function()
  plot(network2, layout=as.matrix(tk_coords),
       vertex.size=20,
       vertex.color=my_color, 
       vertex.label.cex=1.1,
       vertex.label.family="Helvetica",
       vertex.label.color="black",
       vertex.frame.color="transparent"
  )
)

## and a heat

new_x <- x[good_genes,]
va <- columnAnnotation(df = data.frame(Environment=meta$Environment,
                                       Genotype=meta$Genotype),
                       col = list(Genotype = c(AtgW = "darkorange1",
                                               AtgE = "dodgerblue2",
                                               AtgH = "navyblue",
                                               B6 = "darkorange1",
                                               NOD2 = "mediumspringgreen"),
                                  Environment=c(wild = "red3", lab = "mediumorchid3")))
hh=grid.grabExpr(draw(
  Heatmap(t(scale(t(new_x))), top_annotation = va,
          col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
          cluster_columns = T, cluster_rows = T,
          show_row_names = T, show_column_names = F,
          row_names_gp = gpar(fontsize = 10),
          heatmap_legend_param=list(title = ""))
))

pdf("fig4/mod_11.1.2_network.pdf", height = 7, width = 20)
grid.arrange(netdraw, hh, nrow=1)
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
y=as.character(meta$Gender)
y[y=="M"]<-2
y[y=="F"]<-1
y=as.numeric(y)
GSA.obj<-GSA(x,y, genenames = rownames(x), genesets=custom_genesets[[1]], 
             resp.type="Two class unpaired", nperms=100, minsize = 10)
GSA.listsets(GSA.obj, geneset.names=gnames,FDRcut=0.1)
GSA.plot(GSA.obj)

y=as.character(meta$Pregnant)
y[y=="Y"]<-2
y[y=="N"]<-1
y=as.numeric(y)
GSA.obj<-GSA(x,y, genenames = rownames(x), genesets=custom_genesets[[1]], 
             resp.type="Two class unpaired", nperms=100, minsize = 10)
GSA.listsets(GSA.obj, geneset.names=gnames,FDRcut=0.1)
GSA.plot(GSA.obj)


##
## genotype "anova-like" test for any difference
y=as.character(meta$Genotype)
y[y=="AtgH"]<-1
y[y=="AtgE"]<-2
y[y=="B6"]<-3
y[y=="NOD2"]<-4

y=as.numeric(y)
GSA.obj<-GSA(x,y, genenames = rownames(x), genesets=custom_genesets[[1]], 
             resp.type="Multiclass", nperms=100, minsize = 10)
good_gs=rbind(GSA.listsets(GSA.obj, geneset.names=gnames,FDRcut=0.7)$positive, 
              GSA.listsets(GSA.obj, geneset.names=gnames,FDRcut=0.7)$negative)


hey = GSA_networker(good_gs[1:20,1], 
                    x_data = x, 
                    y_data = meta$Genotype, meta_entry=meta,
                    box_cat = "Genotype", 
                    box_cols = c("AtgE"="dodgerblue2", "AtgH"="navyblue",
                                 "B6"="darkorange1", "NOD2"="mediumspringgreen"))

namer <- paste0("genotype/ANOVA_like_comp_nets.pdf")
pdf(namer, height = 25, width = 30)
print(marrangeGrob(grobs=hey,ncol=1,nrow=1))
dev.off()
#
#
#
#
#
#
#
#
### geno GSA split lab and wild

y=as.character(meta$Genotype)
y[y=="AtgH"]<-1
y[y=="AtgE"]<-2
y[y=="B6"]<-3
y[y=="NOD2"]<-4

y=as.numeric(y)

env <- c("lab", "wild")
for(b in 1:length(env)){
  meta_sub2 <- subset(meta, Environment == env[b])
  
  sub_char <- as.character(meta_sub2$mouse_id)
  new_x <- x[,sub_char]
  
  new_y=as.character(meta_sub2$Genotype)
  new_y=as.numeric(as.factor(new_y))
  
  try({
    GSA.obj<-GSA(new_x,new_y, genenames = rownames(new_x), 
                 genesets=custom_genesets[[1]], 
                 resp.type="Multiclass", nperms=100)
    geno_add=as.data.frame(rbind(GSA.listsets(GSA.obj, geneset.names=gnames,FDRcut=0.8)$positive, 
                                 GSA.listsets(GSA.obj, geneset.names=gnames,FDRcut=0.8)$negative))
  })

  hey = GSA_networker(geno_add[,1], 
                      x_data = x[,as.character(meta_sub2$mouse_id)], 
                      y_data = meta_sub2$Genotype, meta_entry=meta_sub2,
                      box_cat = "Genotype", 
                      box_cols = c("AtgE"="dodgerblue2", "AtgH"="navyblue",
                                   "B6"="darkorange1", "NOD2"="mediumspringgreen"))
  
  namer <- paste0("genotype/ANOVA_like_comp", env[b], "_nets.pdf")
  pdf(namer, height = 25, width = 30)
  print(marrangeGrob(grobs=hey,ncol=1,nrow=1))
  dev.off()
}


#
##

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
#
#
#
#
#
#

######### part 2 investigate the important parts of our network!

### IL-5 bacteroides
custom_genesets[[1]][134]

hey = GSA_networker(c(134), 
                    x_data = x, 
                    y_data = meta$Genotype, meta_entry=meta,
                    box_cat = "Genotype", 
                    box_cols = c("AtgE"="dodgerblue2", "AtgH"="navyblue",
                                 "B6"="darkorange1", "NOD2"="mediumspringgreen"))

namer <- paste0("fig4/IL_5_bacteroides_genotype_nets.pdf")
pdf(namer, height = 25, width = 30)
print(marrangeGrob(grobs=hey,ncol=1,nrow=1))
dev.off()

hey = GSA_networker(c(134), 
                    x_data = x, 
                    y_data = meta$Environment, meta_entry=meta,
                    box_cat = "Environment", 
                    box_cols = c("lab"="mediumorchid3", "wild"="red3"))

namer <- paste0("fig4/IL_5_bacteroides_env_nets.pdf")
pdf(namer, height = 25, width = 30)
print(marrangeGrob(grobs=hey,ncol=1,nrow=1))
dev.off()


### what genes are most associated with IL-5

il5_df <- data.frame(net_mat[,"IL.5_BacteroidesV"], names=rownames(net_mat))

# most significant association is 2.2.2
# what genes are in there?

module_key <- read.table(
  "/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/data/MLN_RNA_Seq/gene_module_key.txt", T)

#write a function that plots genes from a gene module according to lab/wild and genotype
#load the real genes

gene_mat <- readRDS("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/IntegrationModel/ml_inputs/genes.RData")
gene_mat2 <- t(gene_mat)
gene_mat2 <- gene_mat2[which(rownames(gene_mat2) %in% as.character(module_key$Genes)),]

module_of_int = "2.1.1.1.1.1.2"

gene_heater <- function(genes_of_int, meta_entry, x_data, x_of_int){
  #genes_of_int <- as.character(subset(module_key, module %in% module_of_int)$Genes)
  new_x <- x_data[genes_of_int,]
  
  x_of_int <- as.matrix(x_of_int)
  corr_corr <- corr.test(log2(t(new_x)+1), x_of_int)
  corr_df <- data.frame(R=corr_corr$r, pval=corr_corr$p, names=rownames(corr_corr$r))
  
  corr_df2 <- subset(corr_df, abs(R) > 0.25)
  corr_df2 <- corr_df2[order(corr_df2$R, decreasing=T),]
  corr_df2$names <- factor(corr_df2$names, levels = as.character(corr_df2$names))
  corr_df2$col <- "down"
  corr_df2$col[corr_df2$R > 0] <- "up"
  
  corr_bar <- ggplot(corr_df2, aes(names, R, fill=col)) +
    geom_col() + coord_flip() +
    scale_fill_manual(values = c(up="red", down="dodgerblue")) +
    theme_bw() + 
    theme(legend.position = 'none',
          axis.text = element_text(color="black"))
  
  gene_filt <- as.character(corr_df2$names)
  ##heat
  new_x <- x_data[rev(gene_filt),]
  va <- columnAnnotation(df = data.frame(Environment=meta_entry$Environment,
                                         Genotype=meta_entry$Genotype),
                         col = list(Genotype = c(AtgW = "darkorange1",
                                                 AtgE = "dodgerblue2",
                                                 AtgH = "navyblue",
                                                 B6 = "darkorange1",
                                                 NOD2 = "mediumspringgreen"),
                                    Environment=c(wild = "red3", lab = "mediumorchid3")))
  hh=grid.grabExpr(draw(
    Heatmap(t(scale(t(new_x))), top_annotation = va,
            col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
            cluster_columns = T, cluster_rows = F,
            show_row_names = T, show_column_names = F,
            row_names_gp = gpar(fontsize = 10),
            heatmap_legend_param=list(title = ""))
  ))
  
  return(arrangeGrob(hh, corr_bar, nrow=1, widths = c(2,1)))
  
}

x_of_int <- unlist(c(x["IL.5_BacteroidesV",]))

genes_of_int <- as.character(subset(module_key, module == "11.1.2" | 
                                      module == "11.1.1.2" | 
                                      module == "11.1.1.1")$Genes)

il5_plots <- gene_heater(genes_of_int, meta, gene_mat2, x_of_int)
pdf("fig4/IL5_bacteroides_comp_to_11.1s.pdf", height = 5.5, width = 10)
grid.arrange(il5_plots)
dev.off()

#
#il5 box plot

il5_box_df <- data.frame(il5=x_of_int, meta)
il5_box = ggplot(il5_box_df, aes(Environment, il5, color= Environment, fill= Environment))+
  geom_boxplot(alpha=0.5, outlier.shape=NA) +
  #geom_jitter(aes(shape=Genotype), width = 0.2) +
  geom_jitter(width = 0.2) +
  coord_cartesian(ylim=c(0,14))+
  ggtitle("IL-5 Bacteroides") +
  scale_color_manual(values = c(wild = "red3", lab = "mediumorchid3")) +
  scale_fill_manual(values = c(wild = "red3", lab = "mediumorchid3")) +
  theme_bw() +
  ylab("Log2 Value") + #xlab(box_cat) +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12, color="black"),
        legend.position = 'none',
        legend.title = element_text(size=15),
        legend.text = element_text(size=12, color="black")
  )
  

wilcox.test(il5 ~ Environment, data=il5_box_df)$p.value

### 11.1 gene modules

mod_11.1_df <- data.frame(t(x["11.1.2",]), t(x["11.1.1.2",]), t(x["11.1.1.1",]))
mod_11.1_df$avg_11.1 <- rowMeans(mod_11.1_df)
mod_11.1_df <- cbind(mod_11.1_df, meta)

mod_11.1_box = ggplot(mod_11.1_df, aes(Environment, log2(avg_11.1+1), 
                                       color= Environment, fill= Environment))+
  geom_boxplot(alpha=0.5, outlier.shape=NA) +
  #geom_jitter(aes(shape=Genotype), width = 0.2) +
  geom_jitter(width = 0.2) +
  coord_cartesian(ylim=c(3.75,5.25))+
  ggtitle("Modules 11.1") +
  scale_color_manual(values = c(wild = "red3", lab = "mediumorchid3")) +
  scale_fill_manual(values = c(wild = "red3", lab = "mediumorchid3")) +
  theme_bw() +
  ylab("Log2 Value") + #xlab(box_cat) +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12, color="black"),
        legend.position = 'none',
        legend.title = element_text(size=15),
        legend.text = element_text(size=12, color="black")
  )

wilcox.test(log2(avg_11.1+1) ~ Environment, data=mod_11.1_df)$p.value

#######
#
#
#####
# make heat of module 11.1

mod_11.1_df <- data.frame(t(x["11.1.2",]), t(x["11.1.1.2",]), t(x["11.1.1.1",]))
mod_11.1_df$avg_11.1 <- rowMeans(mod_11.1_df)
mod_11.1_df <- cbind(mod_11.1_df, meta)

module_key <- read.table(
  "/Volumes/research/lokep01labspace/Rewilding_Data/int/data/MLN_RNA_Seq/gene_module_key.txt", T)

#write a function that plots genes from a gene module according to lab/wild and genotype
#load the real genes

gene_mat <- readRDS("/Volumes/research/lokep01labspace/Rewilding_Data/int/IntegrationModel/ml_inputs/genes.RData")
gene_mat2 <- t(gene_mat)
gene_mat2 <- gene_mat2[which(rownames(gene_mat2) %in% as.character(module_key$Genes)),]

genes_of_int <- as.character(subset(module_key, module == "11.1.2" | 
                                      module == "11.1.1.2" | 
                                      module == "11.1.1.1")$Genes)


new_x <- gene_mat2[genes_of_int,]
va <- columnAnnotation(df = data.frame(Environment=meta$Environment,
                                       Genotype=meta$Genotype),
                       col = list(Genotype = c(AtgW = "darkorange1",
                                               AtgE = "dodgerblue2",
                                               AtgH = "navyblue",
                                               B6 = "darkorange1",
                                               NOD2 = "mediumspringgreen"),
                                  Environment=c(wild = "red3", lab = "mediumorchid3")))
hh=grid.grabExpr(draw(
  Heatmap(t(scale(t(new_x))), top_annotation = va,
          col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
          cluster_columns = T, cluster_rows = T,
          show_row_names = T, show_column_names = F,
          row_names_gp = gpar(fontsize = 10),
          heatmap_legend_param=list(title = ""))
))

grid.arrange(hh)


#### #### which of the 11.1 modules differ by lab and wild
de_df <- subset(
  read.table("/Volumes/research/lokep01labspace/Rewilding_Data/int/data/MLN_RNA_Seq/Lab_v_wild_all.txt",
                    header=T, sep='\t'),
  Gene %in% genes_of_int & padj < 0.05 & abs(log2FoldChange) > 0.6)

namer_df <- data.frame(heats=as.character(rownames(new_x)))
namer_df$heats <- as.character(namer_df$heats)
namer_df$new_name <- 0
namer_df$new_name[namer_df$heats %in% as.character(de_df$Gene)] <- 1
namer_df$heats[namer_df$new_name==0]<-""

new_x2=new_x
rownames(new_x2) <- as.character(namer_df$heats)
hh=grid.grabExpr(draw(
  Heatmap(t(scale(t(new_x2))), top_annotation = va,
          col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
          cluster_columns = T, cluster_rows = T,
          show_row_names = T, show_column_names = F,
          row_names_gp = gpar(fontsize = 10),
          heatmap_legend_param=list(title = ""))
))

grid.arrange(hh)

###genes from panther
panth_genes <- unique(c("Snx3", "Cltc", "Arhgap1", "Tmem173", "Tgtp2", "Tgtp1", "Ifi47", 
                        "Gbp2", "Igtp", "Cd24a", "Selplg", "Stk10", "Sema4d", "Itgb2", 
                        "Syk", "Rasgrp1", "Dock2", "Syk", "Myo1f", "Tmem173", "Tgtp2", 
                        "Tgtp1", "Ifi47", "Gbp2", "Igtp", "H2-DMa", "Gata3", "Dock2", 
                        "Syk", "Stat6", "Gata3", "Dock2", "Elf4", "Gpr18", "Rasgrp1", 
                        "Cd24a", "Ets1", "Il2rg", "Dock8", "Dpp4", "Cd47", "Sash3", 
                        "Ctnnb1", "Prkdc", "Dock11","Cd47", "Ets1", "Vamp8", "Cd55", 
                        "Selplg", "Fcrl1" ,"Git1", "Prkdc"))

hh=grid.grabExpr(draw(
  Heatmap(t(scale(t(new_x[panth_genes,]))), top_annotation = va,
          col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
          cluster_columns = T, cluster_rows = T,
          show_row_names = T, show_column_names = F,
          row_names_gp = gpar(fontsize = 10),
          heatmap_legend_param=list(title = ""))
))

pdf("fig4/mod_11.1_heat.pdf", height = 5.5, width = 7)
grid.arrange(hh)
dev.off()


#
#
#load in panther results with top pathways and p-values
panther_res <- read.table("panther_11.1.txt", T, '\t')
panther_small <- panther_res[1:10,]
panther_small$GO.biological.process.complete <- gsub(" \\(G.*", "", panther_small$GO.biological.process.complete)
panther_small$GO.biological.process.complete <- factor(panther_small$GO.biological.process.complete, levels = rev(as.character(panther_small$GO.biological.process.complete)))
panther_small$Client.Text.Box.Input..fold.Enrichment. <- as.numeric(as.character(panther_small$Client.Text.Box.Input..fold.Enrichment.))

fplot=ggplot(panther_small, aes(GO.biological.process.complete, 
                          Client.Text.Box.Input..fold.Enrichment.)) +
geom_col(fill='navy') + coord_flip() +
  ggtitle("Module 11.1 Pathway Enrichment") +
  scale_color_manual(values="navy") +
  theme_bw() +
  ylab("Fold Enrichment") + xlab("Pathway") +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12, color="black"),
        legend.position = 'none',
        legend.title = element_text(size=15),
        legend.text = element_text(size=12, color="black")
  )


###
panther_res <- read.table("panther_11.1.txt", T, '\t')
panther_small <- panther_res[1:10,]
panther_small$GO.biological.process.complete <- gsub(" \\(G.*", "", panther_small$GO.biological.process.complete)
panther_small=panther_small[order(panther_small$Client.Text.Box.Input..raw.P.value., decreasing=F),]
panther_small$GO.biological.process.complete <- factor(panther_small$GO.biological.process.complete, levels = rev(as.character(panther_small$GO.biological.process.complete)))
panther_small$Client.Text.Box.Input..fold.Enrichment. <- as.numeric(as.character(panther_small$Client.Text.Box.Input..fold.Enrichment.))

pplot=ggplot(panther_small, aes(GO.biological.process.complete, 
                          -log10(Client.Text.Box.Input..raw.P.value.))) +
  geom_col(fill='navy') + coord_flip() +
  ggtitle("Module 11.1 Pathway Enrichment") +
  scale_color_manual(values="navy") +
  theme_bw() +
  ylab("-log10(Pvalue)") + xlab("Pathway") +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12, color="black"),
        legend.position = 'none',
        legend.title = element_text(size=15),
        legend.text = element_text(size=12, color="black")
  )

pdf("fig4/panther_pathways.pdf", height = 4, width = 14)
grid.arrange(fplot, pplot, nrow=1)
dev.off()
#
#

### Blood_CD4_T_cellsCD44hi_CD62Lhi

#
blood_cd4_df <- data.frame(t(x["Blood_CD4_T_cells_CD44hi_CD62Lhi",]), meta)
wilcox.test(Blood_CD4_T_cells_CD44hi_CD62Lhi ~ Environment, data=blood_cd4_df)$p.value

blood_cd4_df = ggplot(blood_cd4_df, aes(Environment, Blood_CD4_T_cells_CD44hi_CD62Lhi, color= Environment, fill= Environment))+
  geom_boxplot(alpha=0.5, outlier.shape=NA) +
  #geom_jitter(aes(shape=Genotype), width = 0.2) +
  geom_jitter(width = 0.2) +
  coord_cartesian(ylim=c(0,12))+
  ggtitle("Blood_CD4_T_cells_CD44hi_CD62Lhi") +
  scale_color_manual(values = c(wild = "red3", lab = "mediumorchid3")) +
  scale_fill_manual(values = c(wild = "red3", lab = "mediumorchid3")) +
  theme_bw() +
  ylab("Log2 Value") + #xlab(box_cat) +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12, color="black"),
        legend.position = 'none',
        legend.title = element_text(size=15),
        legend.text = element_text(size=12, color="black")
  )




pdf("fig4/network_modules_of_int_boxes.pdf", height = 5, width = 12)
grid.arrange(il5_box, mod_11.1_box, blood_cd4_df, nrow=1)
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

#
#
###### plot heat of modules 39 and 22

module_key <- read.table(
  "/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/data/MLN_RNA_Seq/gene_module_key.txt", T)

#write a function that plots genes from a gene module according to lab/wild and genotype
#load the real genes

gene_mat <- readRDS("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/IntegrationModel/ml_inputs/genes.RData")
gene_mat2 <- t(gene_mat)
gene_mat2 <- gene_mat2[which(rownames(gene_mat2) %in% as.character(module_key$Genes)),]


genes_of_int <- as.character(subset(module_key, module == "39")$Genes)


new_x <- gene_mat2[genes_of_int,]
va <- columnAnnotation(df = data.frame(Environment=meta$Environment,
                                       Genotype=meta$Genotype),
                       col = list(Genotype = c(AtgW = "darkorange1",
                                               AtgE = "dodgerblue2",
                                               AtgH = "navyblue",
                                               B6 = "darkorange1",
                                               NOD2 = "mediumspringgreen"),
                                  Environment=c(wild = "red3", lab = "mediumorchid3")))
hh39=grid.grabExpr(draw(
  Heatmap(t(scale(t(new_x))), top_annotation = va,
          col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
          cluster_columns = T, cluster_rows = T,
          show_row_names = T, show_column_names = F,
          row_names_gp = gpar(fontsize = 10),
          heatmap_legend_param=list(title = ""))
))


genes_of_int <- as.character(subset(module_key, module == "22")$Genes)
genes_of_int <- genes_of_int[-36]

new_x <- gene_mat2[genes_of_int,]
va <- columnAnnotation(df = data.frame(Environment=meta$Environment,
                                       Genotype=meta$Genotype),
                       col = list(Genotype = c(AtgW = "darkorange1",
                                               AtgE = "dodgerblue2",
                                               AtgH = "navyblue",
                                               B6 = "darkorange1",
                                               NOD2 = "mediumspringgreen"),
                                  Environment=c(wild = "red3", lab = "mediumorchid3")))
hh22=grid.grabExpr(draw(
  Heatmap(t(scale(t(new_x))), top_annotation = va,
          col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
          cluster_columns = T, cluster_rows = T,
          show_row_names = T, show_column_names = F,
          row_names_gp = gpar(fontsize = 10),
          heatmap_legend_param=list(title = ""))
))


pdf("mods_39_and_22_heats.pdf", height = 10, width = 8)
grid.arrange(hh39, hh22, nrow=1)
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
#
#
#
#

##### genos GSA for each comp

genos = as.character(unique(meta$Genotype))
geno_combo <- combn(unique(meta$Genotype),2)
geno_gs <- as.data.frame(t(good_gs[1,]))
geno_gs$comp=NA
for(i in 1:ncol(geno_combo)){
  g_int1 = as.character(geno_combo[1,i])
  g_int2 = as.character(geno_combo[2,i])
  meta_sub <- subset(meta, Genotype == g_int1 | Genotype == g_int2)
  
  sub_char <- as.character(meta_sub$mouse_id)
  new_x <- x[,sub_char]

  new_y=as.character(meta_sub$Genotype)
  new_y=as.numeric(as.factor(new_y))
  
  try({
    GSA.obj<-GSA(new_x,new_y, genenames = rownames(new_x), genesets=custom_genesets[[1]], resp.type="Two class unpaired", nperms=100)
    GSA.listsets(GSA.obj, geneset.names=gnames,FDRcut=0.1)
    geno_add=as.data.frame(rbind(GSA.listsets(GSA.obj, geneset.names=gnames,FDRcut=0.2)$positive, 
                                 GSA.listsets(GSA.obj, geneset.names=gnames,FDRcut=0.2)$negative))
    geno_add$comp <- paste0(g_int1, "_vs_", g_int2)
    geno_gs <- rbind(geno_gs, geno_add)
  })
}
geno_gs<-geno_gs[-1,]

#
#
#
#
#
#

### make subnetworks!
for ( m in 1:length(unique(geno_gs$comp))){
  gsa_nums <- subset(geno_gs, comp == unique(geno_gs$comp)[m])
  
  meta_sub <- subset(meta, Genotype %in% gsub("_.*", "",unique(geno_gs$comp)[m]) |
                           Genotype %in% gsub(".*_", "",unique(geno_gs$comp)[m]))
  
  hey = GSA_networker(gsa_nums[,1], 
                      x_data = x[,as.character(meta_sub$mouse_id)], 
                      y_data = meta_sub$Genotype, meta_entry=meta_sub,
                      box_cat = "Genotype", 
                      box_cols = c("AtgE"="dodgerblue2", "AtgH"="navyblue",
                                   "B6"="darkorange1", "NOD2"="mediumspringgreen"))
  
  namer <- paste0("genotype/", unique(geno_gs$comp)[m], "_nets.pdf")
  pdf(namer, height = 25, width = 30)
  print(marrangeGrob(grobs=hey,ncol=1,nrow=1))
  dev.off()
  
}

#
#
#
#
#
#
##### lab and wild seems to confound split by enviro and run again!
genos = as.character(unique(meta$Genotype))
geno_combo <- combn(unique(meta$Genotype),2)
geno_gs <- as.data.frame(t(good_gs[1,]))
geno_gs$comp=NA
geno_gs$env <- NA
for(i in 1:ncol(geno_combo)){
  g_int1 = as.character(geno_combo[1,i])
  g_int2 = as.character(geno_combo[2,i])
  meta_sub <- subset(meta, Genotype == g_int1 | Genotype == g_int2)
  
  env <- c("lab", "wild")
  for(b in 1:length(env)){
    meta_sub2 <- subset(meta_sub, Environment == env[b])
    
    sub_char <- as.character(meta_sub2$mouse_id)
    new_x <- x[,sub_char]
    
    new_y=as.character(meta_sub2$Genotype)
    new_y=as.numeric(as.factor(new_y))
    
    try({
      GSA.obj<-GSA(new_x,new_y, genenames = rownames(new_x), genesets=custom_genesets[[1]], resp.type="Two class unpaired", nperms=100)
      GSA.listsets(GSA.obj, geneset.names=gnames,FDRcut=0.2)
      geno_add=as.data.frame(rbind(GSA.listsets(GSA.obj, geneset.names=gnames,FDRcut=0.2)$positive, 
                                   GSA.listsets(GSA.obj, geneset.names=gnames,FDRcut=0.2)$negative))
      geno_add$comp <- paste0(g_int1, "_vs_", g_int2)
      geno_add$env <- env[b]
      geno_gs <- rbind(geno_gs, geno_add)
    })
  }
}
geno_gs<-geno_gs[-1,]


for ( m in 1:length(unique(geno_gs$comp))){
  gsa_nums <- subset(geno_gs, comp == unique(geno_gs$comp)[m])
  
  meta_sub <- subset(meta, Genotype %in% gsub("_.*", "",unique(geno_gs$comp)[m]) |
                       Genotype %in% gsub(".*_", "",unique(geno_gs$comp)[m]))
  
  gsa_lab <- subset(gsa_nums, env == "lab")
  meta_lab <- subset(meta_sub, Environment == "lab")
  
  if(nrow(gsa_lab)>0){
  hey = GSA_networker(gsa_lab[,1], 
                      x_data = x[,as.character(meta_lab$mouse_id)], 
                      y_data = meta_lab$Genotype, meta_entry=meta_lab,
                      box_cat = "Genotype", 
                      box_cols = c("AtgE"="dodgerblue2", "AtgH"="navyblue",
                                   "B6"="darkorange1", "NOD2"="mediumspringgreen"))
  
  namer <- paste0("genotype/lab_", unique(geno_gs$comp)[m], "_nets.pdf")
  pdf(namer, height = 25, width = 30)
  print(marrangeGrob(grobs=hey,ncol=1,nrow=1))
  dev.off()
  }
  
  gsa_wild <- subset(gsa_nums, env == "wild")
  meta_wild <- subset(meta_sub, Environment == "wild")
  
  if(nrow(gsa_wild)>0){
  hey = GSA_networker(gsa_wild[,1], 
                      x_data = x[,as.character(meta_wild$mouse_id)], 
                      y_data = meta_wild$Genotype, meta_entry=meta_wild,
                      box_cat = "Genotype", 
                      box_cols = c("AtgE"="dodgerblue2", "AtgH"="navyblue",
                                   "B6"="darkorange1", "NOD2"="mediumspringgreen"))
  
  namer <- paste0("genotype/wild_", unique(geno_gs$comp)[m], "_nets.pdf")
  pdf(namer, height = 25, width = 30)
  print(marrangeGrob(grobs=hey,ncol=1,nrow=1))
  dev.off()
  }
  
}



#
#
#

#
#
#
### what about candida
grep("IFN.y_Candi", custom_genesets[[1]])

hey = GSA_networker(c("116"), 
                    x_data = x[,as.character(meta$mouse_id)], 
                    y_data = meta$Environment, meta_entry=meta,
                    box_cat = "Environment", 
                    box_cols = c("lab"="mediumorchid3", "wild"="red3"))

namer <- paste0("IFN.y_StaphA_nets.pdf")
pdf(namer, height = 25, width = 30)
print(marrangeGrob(grobs=hey,ncol=1,nrow=1))
dev.off()


#
#
##
### what about all candida

### env
ugh=grep("Candida", custom_genesets[[1]])

hey = GSA_networker(ugh, 
                    x_data = x[,as.character(meta$mouse_id)], 
                    y_data = meta$Environment, meta_entry=meta,
                    box_cat = "Environment", 
                    box_cols = c("lab"="mediumorchid3", "wild"="red3"))

namer <- paste0("candida/env_all_candida_nets.pdf")
pdf(namer, height = 25, width = 30)
print(marrangeGrob(grobs=hey,ncol=1,nrow=1))
dev.off()

### genotype
genos = as.character(unique(meta$Genotype))
geno_combo <- combn(unique(meta$Genotype),2)
geno_gs <- data.frame(g1=NA)
for(i in 1:ncol(geno_combo)){
  g_int1 = as.character(geno_combo[1,i])
  g_int2 = as.character(geno_combo[2,i])
    geno_add <- paste0(g_int1, "_vs_", g_int2)
    geno_gs <- rbind(geno_gs, geno_add)
}
geno_gs<-geno_gs[-1,]

for ( m in 1:length(unique(geno_gs))){
  gsa_nums <- ugh
  
  meta_sub <- subset(meta, Genotype %in% gsub("_.*", "",unique(geno_gs)[m]) |
                       Genotype %in% gsub(".*_", "",unique(geno_gs)[m]))
  
  hey = GSA_networker(gsa_nums, 
                      x_data = x[,as.character(meta_sub$mouse_id)], 
                      y_data = meta_sub$Genotype, meta_entry=meta_sub,
                      box_cat = "Genotype", 
                      box_cols = c("AtgE"="dodgerblue2", "AtgH"="navyblue",
                                   "B6"="darkorange1", "NOD2"="mediumspringgreen"))
  
  namer <- paste0("candida/", unique(geno_gs)[m], "_nets.pdf")
  pdf(namer, height = 25, width = 30)
  print(marrangeGrob(grobs=hey,ncol=1,nrow=1))
  dev.off()
  
}



###env and genotype
for ( m in 1:length(unique(geno_gs))){
  #gsa_nums <- subset(geno_gs, comp == unique(geno_gs)[m])
  
  meta_sub <- subset(meta, Genotype %in% gsub("_.*", "",unique(geno_gs)[m]) |
                       Genotype %in% gsub(".*_", "",unique(geno_gs)[m]))
  
  gsa_lab <- ugh
  meta_lab <- subset(meta_sub, Environment == "lab")
  
  #if(nrow(gsa_lab)>0){
    hey = GSA_networker(gsa_lab, 
                        x_data = x[,as.character(meta_lab$mouse_id)], 
                        y_data = meta_lab$Genotype, meta_entry=meta_lab,
                        box_cat = "Genotype", 
                        box_cols = c("AtgE"="dodgerblue2", "AtgH"="navyblue",
                                     "B6"="darkorange1", "NOD2"="mediumspringgreen"))
    
    namer <- paste0("candida/lab_", unique(geno_gs)[m], "_nets.pdf")
    pdf(namer, height = 25, width = 30)
    print(marrangeGrob(grobs=hey,ncol=1,nrow=1))
    dev.off()
 # }
  
  gsa_wild <- ugh
  meta_wild <- subset(meta_sub, Environment == "wild")
  
  #if(nrow(gsa_wild)>0){
    hey = GSA_networker(gsa_wild, 
                        x_data = x[,as.character(meta_wild$mouse_id)], 
                        y_data = meta_wild$Genotype, meta_entry=meta_wild,
                        box_cat = "Genotype", 
                        box_cols = c("AtgE"="dodgerblue2", "AtgH"="navyblue",
                                     "B6"="darkorange1", "NOD2"="mediumspringgreen"))
    
    namer <- paste0("candida/wild_", unique(geno_gs)[m], "_nets.pdf")
    pdf(namer, height = 25, width = 30)
    print(marrangeGrob(grobs=hey,ncol=1,nrow=1))
    dev.off()
  #}
  
}





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
#
#
#

##
#

#

#
#
#### pick my own genesets from aim 1a
int = custom_genesets
which(rownames(net_mat)=="IL.5_BacteroidesV")

good_genes <- unlist(custom_genesets[[1]][134])

new_ids <- rownames(net_mat_og)[which(rownames(net_mat_og)%in%good_genes)]

new_net <- net_mat_og[new_ids,new_ids]
up_feat <- mats[[2]][new_ids,]

network=graph_from_adjacency_matrix(as.matrix(new_net), weighted=T, mode="undirected", diag=F)

icoord = layout_with_fr(network)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

E(network)$color[E(network)$weight > 0] <- add.alpha('red', alpha=0.2)
E(network)$color[E(network)$weight < 0] <- add.alpha('blue', alpha=0.2)

coul <- c(cytokines="darkseagreen2", genes="lightcoral", cells="dodgerblue2", otus="mediumorchid3")
my_color=coul[up_feat$dataType]
# plot

curr_coord = icoord
new_coord = data.frame(curr_coord, up_feat)
colnames(new_coord) <- c("X", "Y", 'Feature', "DataType")

#library(ggplot2)
#ggplot(new_coord, aes(X,Y,color=DataType)) + 
#  geom_point()

gs = as.character(unique(new_coord$DataType))
gs=c("genes", "cells", "cytokines", "otus")
better_coords <- data.frame(X=NA,Y=NA,Feature=NA,DataType=NA,X2=NA,Y2=NA)
low_bounds_y=c(250,800,250,0)
high_bounds_y=c(700,950,700,200)

low_bounds_x=c(0,300,700,300)
high_bounds_x=c(300,700,1000,700)

for(j in 1:length(gs)){
  curr <- subset(new_coord, DataType == gs[j])
  
  x2 = sample(low_bounds_x[j]:high_bounds_x[j], nrow(curr))
  y2 = sample(low_bounds_y[j]:high_bounds_y[j], nrow(curr), replace = T)
  
  curr$X2=x2
  curr$Y2=y2
  
  better_coords = rbind(better_coords, curr)
}
better_coords=better_coords[-1,]

#ggplot(better_coords, aes(X2,Y2,color=DataType)) + 
#  geom_point()
better_coords=better_coords[rownames(new_coord),]

namer=paste0("il5_bacteroides_network.pdf")
pdf(namer, height = 7, width = 7)
par(bg="white", mar=c(0,0,0,0))
set.seed(4)
plot(network, layout=as.matrix(better_coords[,5:6]),
     vertex.size=10,
     vertex.color=my_color, 
     vertex.label.cex=1.2,
     vertex.label.family="Helvetica",
     vertex.label.color="black",
     vertex.frame.color="transparent"
)
dev.off()

new_x <- x[good_genes,]

namer=paste0("il5_bacteroides_heat.pdf")
pdf(namer, height = 8, width = 15)
va <- columnAnnotation(df = data.frame(Environment=meta$Environment,
                                       Genotype=meta$Genotype),
                       col = list(Genotype = c(AtgW = "darkorange1",
                                               AtgE = "dodgerblue2",
                                               AtgH = "navyblue",
                                               B6 = "darkorange1",
                                               NOD2 = "mediumspringgreen"),
                                  Environment=c(wild = "red3", lab = "mediumorchid3")))
hh=Heatmap(t(scale(t(new_x))), top_annotation = va,
           col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
           cluster_columns = T, cluster_rows = T,
           show_row_names = T, show_column_names = F,
           row_names_gp = gpar(fontsize = 10),
           heatmap_legend_param=list(title = ""))
print(hh)
dev.off()


#
#
#pick more
which(rownames(net_mat)=="IL.10_PseudomonasA")

good_genes <- unlist(custom_genesets[[1]][148])

new_ids <- rownames(net_mat_og)[which(rownames(net_mat_og)%in%good_genes)]

new_net <- net_mat_og[new_ids,new_ids]
up_feat <- mats[[2]][new_ids,]

network=graph_from_adjacency_matrix(as.matrix(new_net), weighted=T, mode="undirected", diag=F)
icoord = layout_with_fr(network)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

E(network)$color[E(network)$weight > 0] <- add.alpha('red', alpha=0.2)
E(network)$color[E(network)$weight < 0] <- add.alpha('blue', alpha=0.2)

coul <- c(cytokines="darkseagreen2", genes="lightcoral", cells="dodgerblue2", otus="mediumorchid3")
my_color=coul[up_feat$dataType]
# plot

curr_coord = icoord
new_coord = data.frame(curr_coord, up_feat)
colnames(new_coord) <- c("X", "Y", 'Feature', "DataType")

#library(ggplot2)
#ggplot(new_coord, aes(X,Y,color=DataType)) + 
#  geom_point()

gs = as.character(unique(new_coord$DataType))
gs=c("genes", "cells", "cytokines", "otus")
better_coords <- data.frame(X=NA,Y=NA,Feature=NA,DataType=NA,X2=NA,Y2=NA)
low_bounds_y=c(250,800,250,0)
high_bounds_y=c(700,950,700,200)

low_bounds_x=c(0,300,700,300)
high_bounds_x=c(300,700,1000,700)

for(j in 1:length(gs)){
  curr <- subset(new_coord, DataType == gs[j])
  
  x2 = sample(low_bounds_x[j]:high_bounds_x[j], nrow(curr))
  y2 = sample(low_bounds_y[j]:high_bounds_y[j], nrow(curr), replace = T)
  
  curr$X2=x2
  curr$Y2=y2
  
  better_coords = rbind(better_coords, curr)
}
better_coords=better_coords[-1,]

ggplot(better_coords, aes(X2,Y2,color=DataType)) + 
  geom_point()
better_coords=better_coords[rownames(new_coord),]

namer=paste0("MLN_CD4_t_cells_network.pdf")
pdf(namer, height = 10, width = 10)
par(bg="white", mar=c(0,0,0,0))
set.seed(4)
plot(network, layout=as.matrix(better_coords[,5:6]),
     vertex.size=10,
     vertex.color=my_color, 
     vertex.label.cex=1.2,
     vertex.label.family="Helvetica",
     vertex.label.color="black",
     vertex.frame.color="transparent"
)
dev.off()

new_x <- x[good_genes,]

namer=paste0("MLN_CD4_t_cells_heat.pdf")
pdf(namer, height = 8, width = 15)
va <- columnAnnotation(df = data.frame(Environment=meta$Environment,
                                       Genotype=meta$Genotype),
                       col = list(Genotype = c(AtgW = "darkorange1",
                                               AtgE = "dodgerblue2",
                                               AtgH = "navyblue",
                                               B6 = "darkorange1",
                                               NOD2 = "mediumspringgreen"),
                                  Environment=c(wild = "red3", lab = "mediumorchid3")))
hh=Heatmap(t(scale(t(new_x))), top_annotation = va,
           col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
           cluster_columns = T, cluster_rows = T,
           show_row_names = T, show_column_names = F,
           row_names_gp = gpar(fontsize = 10),
           heatmap_legend_param=list(title = ""))
print(hh)
dev.off()





#
#
#
#
#
#most connected node
#pick more
which(rownames(net_mat)=="IL.5_CandidaA")

good_genes <- unlist(custom_genesets[[1]][241])

new_ids <- rownames(net_mat_og)[which(rownames(net_mat_og)%in%good_genes)]

new_net <- net_mat_og[new_ids,new_ids]
up_feat <- mats[[2]][new_ids,]

network=graph_from_adjacency_matrix(as.matrix(new_net), weighted=T, mode="undirected", diag=F)
icoord = layout_with_fr(network)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

E(network)$color[E(network)$weight > 0] <- add.alpha('red', alpha=0.2)
E(network)$color[E(network)$weight < 0] <- add.alpha('blue', alpha=0.2)

coul <- c(cytokines="darkseagreen2", genes="lightcoral", cells="dodgerblue2", otus="mediumorchid3")
my_color=coul[up_feat$dataType]
# plot

curr_coord = icoord
new_coord = data.frame(curr_coord, up_feat)
colnames(new_coord) <- c("X", "Y", 'Feature', "DataType")

#library(ggplot2)
#ggplot(new_coord, aes(X,Y,color=DataType)) + 
#  geom_point()

gs = as.character(unique(new_coord$DataType))
gs=c("genes", "cells", "cytokines", "otus")
better_coords <- data.frame(X=NA,Y=NA,Feature=NA,DataType=NA,X2=NA,Y2=NA)
low_bounds_y=c(250,800,250,0)
high_bounds_y=c(700,950,700,200)

low_bounds_x=c(0,300,700,300)
high_bounds_x=c(300,700,1000,700)

for(j in 1:length(gs)){
  curr <- subset(new_coord, DataType == gs[j])
  
  x2 = sample(low_bounds_x[j]:high_bounds_x[j], nrow(curr))
  y2 = sample(low_bounds_y[j]:high_bounds_y[j], nrow(curr), replace = T)
  
  curr$X2=x2
  curr$Y2=y2
  
  better_coords = rbind(better_coords, curr)
}
better_coords=better_coords[-1,]

ggplot(better_coords, aes(X2,Y2,color=DataType)) + 
  geom_point()
better_coords=better_coords[rownames(new_coord),]

namer=paste0("il5_candida_network.pdf")
pdf(namer, height = 10, width = 10)
par(bg="white", mar=c(0,0,0,0))
set.seed(4)
plot(network, layout=as.matrix(better_coords[,5:6]),
     vertex.size=10,
     vertex.color=my_color, 
     vertex.label.cex=1.2,
     vertex.label.family="Helvetica",
     vertex.label.color="black",
     vertex.frame.color="transparent"
)
dev.off()

new_x <- x[good_genes,]

namer=paste0("il5_candida_heat.pdf")
pdf(namer, height = 8, width = 15)
va <- columnAnnotation(df = data.frame(Environment=meta$Environment,
                                       Genotype=meta$Genotype),
                       col = list(Genotype = c(AtgW = "darkorange1",
                                               AtgE = "dodgerblue2",
                                               AtgH = "navyblue",
                                               B6 = "darkorange1",
                                               NOD2 = "mediumspringgreen"),
                                  Environment=c(wild = "red3", lab = "mediumorchid3")))
hh=Heatmap(t(scale(t(new_x))), top_annotation = va,
           col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
           cluster_columns = T, cluster_rows = T,
           show_row_names = T, show_column_names = F,
           row_names_gp = gpar(fontsize = 10),
           heatmap_legend_param=list(title = ""))
print(hh)
dev.off()



#
#
#
#
###
#pick more
which(rownames(net_mat)=="IL.5_BacteroidesV")

good_genes <- unlist(custom_genesets[[1]][240])

new_ids <- rownames(net_mat_og)[which(rownames(net_mat_og)%in%good_genes)]

new_net <- net_mat_og[new_ids,new_ids]
up_feat <- mats[[2]][new_ids,]

network=graph_from_adjacency_matrix(as.matrix(new_net), weighted=T, mode="undirected", diag=F)
icoord = layout_with_fr(network)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

E(network)$color[E(network)$weight > 0] <- add.alpha('red', alpha=0.2)
E(network)$color[E(network)$weight < 0] <- add.alpha('blue', alpha=0.2)

coul <- c(cytokines="darkseagreen2", genes="lightcoral", cells="dodgerblue2", otus="mediumorchid3")
my_color=coul[up_feat$dataType]
# plot

curr_coord = icoord
new_coord = data.frame(curr_coord, up_feat)
colnames(new_coord) <- c("X", "Y", 'Feature', "DataType")

#library(ggplot2)
#ggplot(new_coord, aes(X,Y,color=DataType)) + 
#  geom_point()

gs = as.character(unique(new_coord$DataType))
gs=c("genes", "cells", "cytokines", "otus")
better_coords <- data.frame(X=NA,Y=NA,Feature=NA,DataType=NA,X2=NA,Y2=NA)
low_bounds_y=c(250,800,250,0)
high_bounds_y=c(700,950,700,200)

low_bounds_x=c(0,300,700,300)
high_bounds_x=c(300,700,1000,700)

for(j in 1:length(gs)){
  curr <- subset(new_coord, DataType == gs[j])
  
  x2 = sample(low_bounds_x[j]:high_bounds_x[j], nrow(curr))
  y2 = sample(low_bounds_y[j]:high_bounds_y[j], nrow(curr), replace = T)
  
  curr$X2=x2
  curr$Y2=y2
  
  better_coords = rbind(better_coords, curr)
}
better_coords=better_coords[-1,]

ggplot(better_coords, aes(X2,Y2,color=DataType)) + 
  geom_point()
better_coords=better_coords[rownames(new_coord),]

namer=paste0("il5_bacteroides_network.pdf")
pdf(namer, height = 10, width = 10)
par(bg="white", mar=c(0,0,0,0))
set.seed(4)
plot(network, layout=as.matrix(better_coords[,5:6]),
     vertex.size=10,
     vertex.color=my_color, 
     vertex.label.cex=1.2,
     vertex.label.family="Helvetica",
     vertex.label.color="black",
     vertex.frame.color="transparent"
)
dev.off()

new_x <- x[good_genes,]

namer=paste0("il5_bacteroides_heat.pdf")
pdf(namer, height = 8, width = 15)
va <- columnAnnotation(df = data.frame(Environment=meta$Environment,
                                       Genotype=meta$Genotype),
                       col = list(Genotype = c(AtgW = "darkorange1",
                                               AtgE = "dodgerblue2",
                                               AtgH = "navyblue",
                                               B6 = "darkorange1",
                                               NOD2 = "mediumspringgreen"),
                                  Environment=c(wild = "red3", lab = "mediumorchid3")))
hh=Heatmap(t(scale(t(new_x))), top_annotation = va,
           col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
           cluster_columns = T, cluster_rows = T,
           show_row_names = T, show_column_names = F,
           row_names_gp = gpar(fontsize = 10),
           heatmap_legend_param=list(title = ""))
print(hh)
dev.off()

#
#
#

#gene network loop

ints = c("11.1.1.1", "11.1.2", "22", "1.1.2.1.1.2")
for(k in 1:length(ints)){
  

nummer = which(rownames(net_mat)==ints[k])

good_genes <- unlist(custom_genesets[[1]][nummer])

new_ids <- rownames(net_mat_og)[which(rownames(net_mat_og)%in%good_genes)]

new_net <- net_mat_og[new_ids,new_ids]
up_feat <- mats[[2]][new_ids,]

network=graph_from_adjacency_matrix(as.matrix(new_net), weighted=T, mode="undirected", diag=F)
icoord = layout_with_fr(network)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

E(network)$color[E(network)$weight > 0] <- add.alpha('red', alpha=0.2)
E(network)$color[E(network)$weight < 0] <- add.alpha('blue', alpha=0.2)

coul <- c(cytokines="darkseagreen2", genes="lightcoral", cells="dodgerblue2", otus="mediumorchid3")
my_color=coul[up_feat$dataType]
# plot

curr_coord = icoord
new_coord = data.frame(curr_coord, up_feat)
colnames(new_coord) <- c("X", "Y", 'Feature', "DataType")

#library(ggplot2)
#ggplot(new_coord, aes(X,Y,color=DataType)) + 
#  geom_point()

gs = as.character(unique(new_coord$DataType))
gs=c("genes", "cells", "cytokines", "otus")
better_coords <- data.frame(X=NA,Y=NA,Feature=NA,DataType=NA,X2=NA,Y2=NA)
low_bounds_y=c(250,800,250,0)
high_bounds_y=c(700,950,700,200)

low_bounds_x=c(0,300,700,300)
high_bounds_x=c(300,700,1000,700)

for(j in 1:length(gs)){
  curr <- subset(new_coord, DataType == gs[j])
  
  x2 = sample(low_bounds_x[j]:high_bounds_x[j], nrow(curr))
  y2 = sample(low_bounds_y[j]:high_bounds_y[j], nrow(curr), replace = T)
  
  curr$X2=x2
  curr$Y2=y2
  
  better_coords = rbind(better_coords, curr)
}
better_coords=better_coords[-1,]

ggplot(better_coords, aes(X2,Y2,color=DataType)) + 
  geom_point()
better_coords=better_coords[rownames(new_coord),]

namer=paste0(ints[k], "_network.pdf")
pdf(namer, height = 10, width = 10)
par(bg="white", mar=c(0,0,0,0))
set.seed(4)
plot(network, layout=as.matrix(better_coords[,5:6]),
     vertex.size=10,
     vertex.color=my_color, 
     vertex.label.cex=1.2,
     vertex.label.family="Helvetica",
     vertex.label.color="black",
     vertex.frame.color="transparent"
)
dev.off()

new_x <- x[good_genes,]

namer=paste0(ints[k], "_heat.pdf")
pdf(namer, height = 8, width = 15)
va <- columnAnnotation(df = data.frame(Environment=meta$Environment,
                                       Genotype=meta$Genotype),
                       col = list(Genotype = c(AtgW = "darkorange1",
                                               AtgE = "dodgerblue2",
                                               AtgH = "navyblue",
                                               B6 = "darkorange1",
                                               NOD2 = "mediumspringgreen"),
                                  Environment=c(wild = "red3", lab = "mediumorchid3")))
hh=Heatmap(t(scale(t(new_x))), top_annotation = va,
           col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
           cluster_columns = T, cluster_rows = T,
           show_row_names = T, show_column_names = F,
           row_names_gp = gpar(fontsize = 10),
           heatmap_legend_param=list(title = ""))
print(hh)
dev.off()
}
#
#
#
#
#
#
#
#
### make subnetworks!
geno_gs
sampler = sample(1:nrow(geno_gs), 10)
#for(i in 1:nrow(geno_gs)){
for(i in sampler){
  good_genes <- custom_genesets[[1]][[as.numeric(as.character(geno_gs[i,1]))]]
  
  if(length(good_genes)> 1){
  new_ids <- rownames(net_mat_og)[which(rownames(net_mat_og)%in%good_genes)]
  
  new_net <- as.matrix(net_mat_og[new_ids,new_ids])
  up_feat <- mats[[2]][new_ids,]
  
  network=graph_from_adjacency_matrix(new_net, weighted=T, mode="undirected", diag=T)
  
  icoord = layout_with_fr(network)
  
  add.alpha <- function(col, alpha=1){
    if(missing(col))
      stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2, 
          function(x) 
            rgb(x[1], x[2], x[3], alpha=alpha))  
  }
  
  E(network)$color[E(network)$weight > 0] <- add.alpha('red', alpha=0.2)
  E(network)$color[E(network)$weight < 0] <- add.alpha('blue', alpha=0.2)
  
  #up_feat <- feats_df[rownames(new_net),]
  my_color=coul[as.numeric(as.factor(up_feat$dataType))]
  # plot
  
  curr_coord = icoord
  new_coord = data.frame(curr_coord, up_feat)
  colnames(new_coord) <- c("X", "Y", 'Feature', "DataType")
  
  #library(ggplot2)
  #ggplot(new_coord, aes(X,Y,color=DataType)) + 
  #  geom_point()
  
  gs = as.character(unique(new_coord$DataType))
  gs=c("genes", "cells", "cytokines", "otus")
  better_coords <- data.frame(X=NA,Y=NA,Feature=NA,DataType=NA,X2=NA,Y2=NA)
  low_bounds_y=c(250,800,250,0)
  high_bounds_y=c(700,950,700,200)
  
  low_bounds_x=c(0,300,700,300)
  high_bounds_x=c(300,700,1000,700)
  
  for(j in 1:length(gs)){
    curr <- subset(new_coord, DataType == gs[j])
    
    x2 = sample(low_bounds_x[j]:high_bounds_x[j], nrow(curr))
    y2 = sample(low_bounds_y[j]:high_bounds_y[j], nrow(curr), replace = T)
    
    curr$X2=x2
    curr$Y2=y2
    
    better_coords = rbind(better_coords, curr)
  }
  better_coords=better_coords[-1,]
  
  #ggplot(better_coords, aes(X2,Y2,color=DataType)) + 
  #  geom_point()
  better_coords=better_coords[rownames(new_coord),]
  
  namer=paste0("genotype/top_network_num_", i, ".pdf")
  pdf(namer, height = 20, width = 20)
  par(bg="white", mar=c(0,0,0,0))
  set.seed(4)
  plot(network, layout=as.matrix(better_coords[,5:6]),
       vertex.size=10,
       vertex.color=my_color, 
       vertex.label.cex=1.2,
       vertex.label.family="Helvetica",
       vertex.label.color="black",
       vertex.frame.color="transparent"
  )
  dev.off()
  
  new_x <- x[good_genes,]
  
  
  namer=paste0("genotype/top_network_heat_num_", i, ".pdf")
  pdf(namer, height = 8, width = 15)
  va <- columnAnnotation(df = data.frame(Environment=meta$Environment,
                         Genotype=meta$Genotype),
                         col = list(Genotype = c(AtgW = "darkorange1",
                                                 AtgE = "dodgerblue2",
                                                 AtgH = "navyblue",
                                                 B6 = "darkorange1",
                                                 NOD2 = "mediumspringgreen"),
                                    Environment=c(wild = "red3", lab = "mediumorchid3")))
  hh=Heatmap(t(scale(t(new_x))), top_annotation = va,
             col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
             cluster_columns = T, cluster_rows = T,
             show_row_names = T, show_column_names = F,
             row_names_gp = gpar(fontsize = 6),
             heatmap_legend_param=list(title = ""))
  print(hh)
  dev.off()
  }
  
}



#
#
#
#
#
#
#

#
#
### for genotypes split between lab and wild, enviro is too strong of a confounder!

genos = as.character(unique(meta$Genotype))
geno_combo <- combn(unique(meta$Genotype),2)
lab_geno_gs <- as.data.frame(t(good_gs[1,]))
lab_geno_gs$comp=NA
for(i in 1:ncol(geno_combo)){
  g_int1 = as.character(geno_combo[1,i])
  g_int2 = as.character(geno_combo[2,i])
  meta_sub <- subset(meta, Genotype == g_int1 | Genotype == g_int2)
  meta_sub <- subset(meta_sub, Environment == "lab")
  
  sub_char <- as.character(meta_sub$mouse_id)
  new_x <- x[,sub_char]
  new_x = t(scale(t(new_x)))
  new_x[is.na(new_x)] <- 0
  
  new_y=as.character(meta_sub$Genotype)
  new_y=as.numeric(as.factor(new_y))
  
  try({
    GSA.obj<-GSA(new_x,new_y, genenames = rownames(new_x), genesets=custom_genesets[[1]], 
                 resp.type="Two class unpaired", nperms=100, minsize = 10)
    GSA.listsets(GSA.obj, geneset.names=gnames,FDRcut=0.5)
    geno_add=as.data.frame(rbind(GSA.listsets(GSA.obj, geneset.names=gnames,FDRcut=0.5)$positive, 
                                 GSA.listsets(GSA.obj, geneset.names=gnames,FDRcut=0.5)$negative))
    geno_add$comp <- paste0(g_int1, "_vs_", g_int2)
    lab_geno_gs <- rbind(lab_geno_gs, geno_add)
  })
}
lab_geno_gs<-lab_geno_gs[-1,]



### make subnetworks!
lab_geno_gs
sampler = sample(1:nrow(lab_geno_gs), 20)
#for(i in 1:nrow(geno_gs)){
for(i in sampler){
  good_genes <- custom_genesets[[1]][[as.numeric(lab_geno_gs[i,1])]]
  
  if(length(good_genes)> 1){
    new_ids <- rownames(net_mat_og)[which(rownames(net_mat_og)%in%good_genes)]
    
    new_net <- net_mat_og[new_ids,new_ids]
    up_feat <- mats[[2]][new_ids,]
    
    network=graph_from_adjacency_matrix(new_net, weighted=T, mode="undirected", diag=T)
    
    icoord = layout_with_fr(network)
    
    add.alpha <- function(col, alpha=1){
      if(missing(col))
        stop("Please provide a vector of colours.")
      apply(sapply(col, col2rgb)/255, 2, 
            function(x) 
              rgb(x[1], x[2], x[3], alpha=alpha))  
    }
    
    E(network)$color[E(network)$weight > 0] <- add.alpha('red', alpha=0.2)
    E(network)$color[E(network)$weight < 0] <- add.alpha('blue', alpha=0.2)
    
    #up_feat <- feats_df[rownames(new_net),]
    my_color=coul[as.numeric(as.factor(up_feat$dataType))]
    # plot
    
    curr_coord = icoord
    new_coord = data.frame(curr_coord, up_feat)
    colnames(new_coord) <- c("X", "Y", 'Feature', "DataType")
    
    #library(ggplot2)
    #ggplot(new_coord, aes(X,Y,color=DataType)) + 
    #  geom_point()
    
    gs = as.character(unique(new_coord$DataType))
    gs=c("genes", "cells", "cytokines", "otus")
    better_coords <- data.frame(X=NA,Y=NA,Feature=NA,DataType=NA,X2=NA,Y2=NA)
    low_bounds_y=c(250,800,250,0)
    high_bounds_y=c(700,950,700,200)
    
    low_bounds_x=c(0,300,700,300)
    high_bounds_x=c(300,700,1000,700)
    
    for(j in 1:length(gs)){
      curr <- subset(new_coord, DataType == gs[j])
      
      x2 = sample(low_bounds_x[j]:high_bounds_x[j], nrow(curr))
      y2 = sample(low_bounds_y[j]:high_bounds_y[j], nrow(curr), replace = T)
      
      curr$X2=x2
      curr$Y2=y2
      
      better_coords = rbind(better_coords, curr)
    }
    better_coords=better_coords[-1,]
    
    #ggplot(better_coords, aes(X2,Y2,color=DataType)) + 
    #  geom_point()
    better_coords=better_coords[rownames(new_coord),]
    
    namer=paste0("lab_genotype/top_network_num_", i, ".pdf")
    pdf(namer, height = 20, width = 20)
    par(bg="white", mar=c(0,0,0,0))
    set.seed(4)
    plot(network, layout=as.matrix(better_coords[,5:6]),
         vertex.size=10,
         vertex.color=my_color, 
         vertex.label.cex=1.2,
         vertex.label.family="Helvetica",
         vertex.label.color="black",
         vertex.frame.color="transparent"
    )
    dev.off()
    
    lab_only <- subset(meta, Environment=='lab')
    new_x <- x[good_genes,lab_only$mouse_id]
    
    
    namer=paste0("lab_genotype/top_network_heat_num_", i, ".pdf")
    pdf(namer, height = 8, width = 15)
    va <- columnAnnotation(df = data.frame(Environment=lab_only$Environment,
                                           Genotype=lab_only$Genotype),
                           col = list(Genotype = c(AtgW = "darkorange1",
                                                   AtgE = "dodgerblue2",
                                                   AtgH = "navyblue",
                                                   B6 = "darkorange1",
                                                   NOD2 = "mediumspringgreen"),
                                      Environment=c(wild = "red3", lab = "mediumorchid3")))
    hh=Heatmap(t(scale(t(new_x))), top_annotation = va,
               col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
               cluster_columns = T, cluster_rows = T,
               show_row_names = T, show_column_names = F,
               row_names_gp = gpar(fontsize = 6),
               heatmap_legend_param=list(title = ""))
    print(hh)
    dev.off()
  }
  
}


#
#
#

#
#
#
#special pics

net_mat <- mats[[1]]
##filter
net_mat[net_mat<0.4 & net_mat > -0.4]=0
good_cols <- which(rowSums(abs(net_mat))>0)
net_mat=net_mat[good_cols,good_cols]

feats_df <- mats[[2]][good_cols,]

gene_good = "9.2.2"
which(colnames(net_mat)=="9.2.2")
net_gene <- net_mat[287,]

goodies <- as.numeric(which(net_gene!=0))

small_net <- net_mat[,c(287,goodies)]
