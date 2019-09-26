#### network creator

setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/IntegrationModel/network")

library(igraph)
library(ComplexHeatmap)
library(circlize)

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

otu_cols <- grep("D_0__Bacteria", colnames(net_mat))
colnames(net_mat)[otu_cols] <- gsub(".*__", "", colnames(net_mat)[otu_cols])

#make facs names shorter
mln_cols <- grep("CD44", colnames(net_mat))
colnames(net_mat)[mln_cols] <- gsub("T_cells_", "T_cells\n", colnames(net_mat)[mln_cols])

mln_cols <- grep("MHCII", colnames(net_mat))
colnames(net_mat)[mln_cols] <- gsub("CD11c_", "CD11c\n", colnames(net_mat)[mln_cols])

#
#

#



##

network=graph_from_adjacency_matrix(as.matrix(net_mat), weighted=T, mode="undirected", diag=F)
icoord = layout_with_fr(network)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

E(network)$color[E(network)$weight > 0] <- add.alpha('red', alpha=0.3)
E(network)$color[E(network)$weight < 0] <- add.alpha('blue', alpha=0.3)

coul <- c(cytokines="darkseagreen2", genes="lightcoral", cells="dodgerblue2", otus="mediumorchid3")
my_color=coul[feats_df$dataType]
# plot

pdf("int_network5.26.pdf", height = 50, width = 50)
par(bg="white", mar=c(0,0,0,0))
set.seed(4)
plot(network, layout=icoord, 
     vertex.size=10,
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

pdf("int_network5.26.pdf", height = 50, width = 50)
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

tkplot(network, layout=as.matrix(better_coords[,5:6]),
       vertex.size=5,
       vertex.color=my_color, 
       vertex.label.cex=1.2,
       vertex.label.family="Helvetica",
       vertex.label.color="black",
       vertex.frame.color="grey90"
)


tk_coords = tk_coords(3)

pdf("int_network.pdf", height = 50, width = 70)
par(bg="white", mar=c(0,0,0,0))
set.seed(4)
plot(network, #layout=as.matrix(better_coords[,5:6]),
     layout=tk_coords,
     vertex.size=10,
     vertex.color=my_color, 
     vertex.label.cex=1.2,
     vertex.label.family="Helvetica",
     vertex.label.color="black",
     vertex.frame.color="transparent"
)

legend(x=-1.5, y=-0.6, 
       legend=paste( levels(as.factor(feats_df$DataType[good_cols]))), 
       col = coul , bty = "n", pch=20 , pt.cex = 14, cex = 9, text.col="black" , horiz = F)
dev.off()


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
###
# GSA
genesets=list()
for(i in 1:nrow(net_mat)){
  curr = net_mat[i,]
  gs=c(rownames(net_mat)[i],rownames(net_mat)[which(curr != 0)])
  
  if(length(gs)>3){
    print(length(gs))
    genesets[[i]]<-gs
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
GSA.listsets(GSA.obj, geneset.names=gnames,FDRcut=0.7)
GSA.plot(GSA.obj)
good_gs=rbind(GSA.listsets(GSA.obj, geneset.names=gnames,FDRcut=0.7)$positive, 
          GSA.listsets(GSA.obj, geneset.names=gnames,FDRcut=0.7)$negative)

y=as.character(meta$Gender)
y[y=="M"]<-2
y[y=="F"]<-1
y=as.numeric(y)
GSA.obj<-GSA(x,y, genenames = rownames(x), genesets=custom_genesets[[1]], 
             resp.type="Two class unpaired", nperms=100, minsize = 10)
GSA.listsets(GSA.obj, geneset.names=gnames,FDRcut=0.7)
GSA.plot(GSA.obj)

y=as.character(meta$Pregnant)
y[y=="Y"]<-2
y[y=="N"]<-1
y=as.numeric(y)
GSA.obj<-GSA(x,y, genenames = rownames(x), genesets=custom_genesets[[1]], 
             resp.type="Two class unpaired", nperms=100, minsize = 10)
GSA.listsets(GSA.obj, geneset.names=gnames,FDRcut=0.5)
GSA.plot(GSA.obj)

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
    GSA.listsets(GSA.obj, geneset.names=gnames,FDRcut=0.5)
    geno_add=as.data.frame(rbind(GSA.listsets(GSA.obj, geneset.names=gnames,FDRcut=0.5)$positive, 
                                 GSA.listsets(GSA.obj, geneset.names=gnames,FDRcut=0.5)$negative))
    geno_add$comp <- paste0(g_int1, "_vs_", g_int2)
    geno_gs <- rbind(geno_gs, geno_add)
  })
}
geno_gs<-geno_gs[-1,]

### make subnetworks!
good_gs
for(i in 1:nrow(good_gs)){
  good_genes <- custom_genesets[[1]][[as.numeric(good_gs[i,1])]]
  
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
  
  namer=paste0("lab_wild/top_network_num_", i, ".pdf")
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
  
  namer=paste0("lab_wild/top_network_heat_num_", i, ".pdf")
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


#
#
#### pick my own genesets from aim 1a
int = custom_genesets
which(rownames(net_mat)=="9.2.2")

good_genes <- unlist(custom_genesets[[1]][152])

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

namer=paste0("9.2.2_network.pdf")
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

namer=paste0("9.2.2_heat.pdf")
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
which(rownames(net_mat)=="MLN_CD4_T_cells")

good_genes <- unlist(custom_genesets[[1]][258])

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
sampler = sample(1:nrow(geno_gs), 30)
#for(i in 1:nrow(geno_gs)){
for(i in sampler){
  good_genes <- custom_genesets[[1]][[as.numeric(geno_gs[i,1])]]
  
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
