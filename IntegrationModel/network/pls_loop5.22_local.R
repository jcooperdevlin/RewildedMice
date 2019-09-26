### pls loop

#library(mixOmics)
setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/IntegrationModel/network")

library(compositions)

genes = readRDS("../ml_inputs/gene_modules.RData")
cytokines = readRDS("../ml_inputs/cytokines.RData")
cells = readRDS("../ml_inputs/facs.RData")
#otus = readRDS("../ml_inputs/otus_relab.RData")
otus = readRDS("../ml_inputs/otus.RData")

library(ComplexHeatmap)
library(circlize)

input_mat <- cbind(log2(genes+1),log2(cytokines+1),log2(cells+1),log2(otus+1))

png("input_log2_heat.png", height=20, width = 25, units='in', res=300)
Heatmap(input_mat,
        col = colorRamp2(c(0, 10), c("white", "purple"))
)
dev.off()


input_mat <- cbind(clr(genes),clr(cytokines),clr(cells),clr(otus))

png("input_clr_heat.png", height=20, width = 25, units='in', res=300)
Heatmap(input_mat,
        col = colorRamp2(c(-7, 0, 7), c("blue", "white", "red"))
        )
dev.off()

inputs = list(genes=clr(genes), cytokines=clr(cytokines), cells=clr(cells), otus=clr(otus))
combos = combn(4,2)

#### RAN ON CLUSTER FOR MIX OMICS PLUS P-VALUE NULL DISTR LOOP
# for ( i in 1:ncol(combos)){
#   X=inputs[[combos[1,i]]]
#   X_name=names(inputs)[combos[1,i]]
#   Y=inputs[[combos[2,i]]]
#   Y_name=names(inputs)[combos[2,i]]
#   
#   first.spls <- spls(X, Y, ncomp = 10)
#   first.val <- perf(first.spls, validation = "Mfold", folds = 5, repeats = 10)
#   ncomp = which(first.val$Q2 < 0.0975)[1]
#   
#   true.spls <- spls(X, Y, ncomp = ncomp)
#   
#   cord.X = cor(true.spls$X, true.spls$variates$X[, 1:ncomp], use = "pairwise")
#   cord.Y = cor(true.spls$Y, true.spls$variates$Y[, 1:ncomp], use = "pairwise")
#   true_sim = cord.X %*% t(cord.Y)
#   
#   true_sim <- data.frame(item1=rownames(true_sim), true_sim)
#   namer <- paste0("sim_mats/", X_name, "_vs_", Y_name, ".txt")
#   write.table(true_sim, namer, sep='\t', row.names=F, quote=F)
#   
#   true_flat <- c(cord.X %*% t(cord.Y))
#   true_flat[is.na(true_flat)] <- 0
#   
#   num_sim <- 10000
#   null_mat <- matrix(NA, 1, ncol(X)*ncol(Y))
#   for(j in 1:num_sim){
#     null.spls <- spls(X, Y[sample(rownames(Y)),], ncomp = ncomp)
#     cord.X = cor(null.spls$X, null.spls$variates$X[, 1:ncomp], use = "pairwise")
#     cord.Y = cor(null.spls$Y, null.spls$variates$Y[, 1:ncomp], use = "pairwise")
#     null_sim = cord.X %*% t(cord.Y)
#     
#     null_flat <- c(null_sim)
#     null_mat <- rbind(null_mat, null_flat)
#   }
#   null_mat <- null_mat[-1,]
#   null_mat[is.na(null_mat)] <- 0
#   
#   p.value = NULL
#   for (k in 1:ncol(null_mat)){
#     p=mean(abs(null_mat[,k])>=abs(true_flat[k]))
#     p.value = c(p.value,p)
#   }
#   p.adj = p.adjust(p.value, "BH")
#   p_mat <- matrix(p.adj, ncol(X), ncol(Y))
#   rownames(p_mat) <- colnames(X)
#   colnames(p_mat) <- colnames(Y)
#   
#   p_mat <- data.frame(item1=rownames(p_mat), p_mat)
#   namer <- paste0("sim_mats/", X_name, "_vs_", Y_name, "_p.adj.txt")
#   write.table(p_mat, namer, sep='\t', row.names=F, quote=F)
# }



#
#
#
#####
## now read in results and determine significant network connections!

p_mats <- list.files("sim_mats", full.names = T, pattern = c("p.adj"))
c_mats <- setdiff(list.files("sim_mats", pattern = c("txt"), full.names = T),p_mats)

source("sim_mat_func.R")

yes = sim_combine(c_mats, p_mats)

new_c_mats <- list()
for(i in 1:length(c_mats)){
  c_mats1 = read.table(c_mats[i], T, "\t")
  rownames(c_mats1) <- c_mats1[,1]
  c_mats1=c_mats1[,-1]
  p_mats1 = read.table(p_mats[i], T, "\t")
  rownames(p_mats1) <- p_mats1[,1]
  p_mats1=p_mats1[,-1]
  
  c_flat = c(unlist(c_mats1))
  p_flat = c(unlist(p_mats1))
  c_flat[p_flat>0.05]<-0
  c_mats2 = matrix(c_flat, nrow(c_mats1), ncol(c_mats1))
  rownames(c_mats2) <- rownames(c_mats1)
  colnames(c_mats2) <- colnames(c_mats1)
  new_c_mats[[i]] <- c_mats1
}

spacer=matrix(0,210,210)
rownames(spacer)=rownames(new_c_mats[[4]])
colnames(spacer)=rownames(new_c_mats[[4]])
gene_by_facs_cyts_otus_genes <- cbind(new_c_mats[[4]], 
                                new_c_mats[[5]], 
                                new_c_mats[[6]],
                                spacer)

spacer=matrix(0,36,36)
rownames(spacer)=rownames(new_c_mats[[1]])
colnames(spacer)=rownames(new_c_mats[[1]])
adder=cbind(spacer,
      t(new_c_mats[[2]]),
      new_c_mats[[1]],
      t(new_c_mats[[4]])
)
colnames(adder) <- gsub("\\+", "\\.", colnames(adder))
gene_facs_by_facs_cyts_otus_genes <- rbind(gene_by_facs_cyts_otus_genes, adder)


spacer=matrix(0,104,104)
rownames(spacer)=rownames(new_c_mats[[2]])
colnames(spacer)=rownames(new_c_mats[[2]])
adder=cbind(new_c_mats[[2]],
            spacer,
            new_c_mats[[3]],
            t(new_c_mats[[5]])
)
colnames(adder) <- gsub("\\+", "\\.", colnames(adder))
gene_facs_cyts_by_facs_cyts_otus_genes <- rbind(gene_facs_by_facs_cyts_otus_genes, adder)


spacer=matrix(0,164,164)
rownames(spacer)=colnames(new_c_mats[[3]])
colnames(spacer)=colnames(new_c_mats[[3]])
adder=cbind(t(new_c_mats[[6]]),
            t(new_c_mats[[1]]),
            t(new_c_mats[[3]]),
            spacer
)
colnames(adder) <- gsub("\\+", "\\.", colnames(adder))
gene_facs_cyts_otus_by_facs_cyts_otus_genes <- rbind(gene_facs_cyts_by_facs_cyts_otus_genes, adder)

library(ComplexHeatmap)
library(circlize)

combined_c_mat <- gene_facs_cyts_otus_by_facs_cyts_otus_genes
combined_c_mat[is.na(combined_c_mat)] <- 0
rower <- rowSums(combined_c_mat)

#png("test_network_heat.png", height=20, width = 25, units='in', res=300)
#Heatmap(combined_c_mat,
#        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
#        )
#dev.off()


#
#
#
##network!!!!

library(igraph)
#library(qgraph)
c_mat_net <- as.matrix(combined_c_mat)
rownames(c_mat_net) = gsub("CD25\\.", "CD25\\+", rownames(c_mat_net))
colnames(c_mat_net) = gsub("CD25\\.", "CD25\\+", colnames(c_mat_net))
rownames(c_mat_net) = gsub("CD45\\.", "CD45\\+", rownames(c_mat_net))
colnames(c_mat_net) = gsub("CD45\\.", "CD45\\+", colnames(c_mat_net))
c_mat_net <- c_mat_net[rownames(c_mat_net), rownames(c_mat_net)]
#fix OTUs
otu_cols <- grep("D_0__Bacteria", colnames(c_mat_net))
colnames(c_mat_net)[otu_cols] <- gsub(".*__", "", colnames(c_mat_net)[otu_cols])

#make facs names shorter
mln_cols <- grep("CD44", colnames(c_mat_net))
colnames(c_mat_net)[mln_cols] <- gsub("T_cells_", "T_cells\n", colnames(c_mat_net)[mln_cols])

mln_cols <- grep("MHCII", colnames(c_mat_net))
colnames(c_mat_net)[mln_cols] <- gsub("CD11c_", "CD11c\n", colnames(c_mat_net)[mln_cols])

#
#

#
coul <- c("darkseagreen2", "lightcoral", "dodgerblue2", "mediumorchid3")
feats_df <- data.frame(Feature=rownames(c_mat_net), 
                       DataType=c(rep("Gene", 210), rep("ImmunePop", 36), rep("Cytokine", 104), rep("OTU", 164)))
rownames(feats_df) <- feats_df$Feature
feats_df=feats_df[rownames(c_mat_net),]
##filter
c_mat_net[c_mat_net<0.55 & c_mat_net > -0.55]=0
good_cols <- which(rowSums(abs(c_mat_net))>0)
c_mat_net=c_mat_net[good_cols,good_cols]

##

network=graph_from_adjacency_matrix(c_mat_net, weighted=T, mode="undirected", diag=T)

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

my_color=coul[as.numeric(as.factor(feats_df$DataType[good_cols]))]
# plot

pdf("int_network.pdf", height = 50, width = 50)
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
new_coord = data.frame(curr_coord, feats_df[good_cols,])
colnames(new_coord) <- c("X", "Y", 'Feature', "DataType")

library(ggplot2)
ggplot(new_coord, aes(X,Y,color=DataType)) + 
  geom_point()

gs = as.character(unique(new_coord$DataType))
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

pdf("int_network.pdf", height = 50, width = 50)
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


tk_coords = tk_coords(8)

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

#
#

#
#
### we have our network now build a gmt list for our subnetwork analysis
genesets=list()
for(i in 1:nrow(c_mat_net)){
  curr = c_mat_net[i,]
  gs=c(rownames(c_mat_net)[i],rownames(c_mat_net)[which(curr != 0)])
  
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

x=data.frame(rbind(t(genes), t(cytokines), t(cells), t(otus)))
meta <- readRDS("../ml_inputs/mouse_metadata.RData")
y=as.character(meta$Environment)
y[y=="lab"]<-2
y[y=="wild"]<-1
y=as.numeric(y)

GSA.obj<-GSA(x,y, genenames = rownames(x), genesets=custom_genesets[[1]], resp.type="Two class unpaired", nperms=100)
GSA.listsets(GSA.obj, geneset.names=geneset.names,FDRcut=0.5)
GSA.plot(GSA.obj)
good_gs=GSA.listsets(GSA.obj, geneset.names=geneset.names,FDRcut=0.5)$positive

y=as.character(meta$Gender)
y[y=="M"]<-2
y[y=="F"]<-1
y=as.numeric(y)
GSA.obj<-GSA(x,y, genenames = rownames(x), genesets=custom_genesets[[1]], resp.type="Two class unpaired", nperms=100)
GSA.listsets(GSA.obj, geneset.names=geneset.names,FDRcut=0.5)
GSA.plot(GSA.obj)


y=as.character(meta$Genotype)
y=as.numeric(as.factor(y))

GSA.obj<-GSA(x,y, genenames = rownames(x), genesets=custom_genesets[[1]], resp.type="Multiclass", nperms=100)
GSA.listsets(GSA.obj, geneset.names=geneset.names,FDRcut=0.5)
GSA.plot(GSA.obj)


### make subnetworks!
good_gs
for(i in 1:20){
  good_genes <- custom_genesets[[1]][[as.numeric(good_gs[i,1])]]
  
  new_ids <- which(rownames(c_mat_net)%in%good_genes)
  
  new_net <- c_mat_net[new_ids,new_ids]
  
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
  
  up_feat <- feats_df[rownames(new_net),]
  my_color=coul[as.numeric(as.factor(up_feat$DataType))]
  # plot
  
  curr_coord = icoord
  new_coord = data.frame(curr_coord, up_feat)
  colnames(new_coord) <- c("X", "Y", 'Feature', "DataType")
  
  #library(ggplot2)
  #ggplot(new_coord, aes(X,Y,color=DataType)) + 
  #  geom_point()
  
  gs = as.character(unique(new_coord$DataType))
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
  
  namer=paste0("top_network_num_", i, ".pdf")
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
  
  
}
