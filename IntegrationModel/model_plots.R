###### Classifier importance and naive networks

### 6.28

library(igraph)
library(ggplot2)
library(gridExtra)
library(grid)
library(ggplotify)
library(reshape)
library(caret)
library(matrixStats)

setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/IntegrationModel")

mod_importance <- function(model, feat_vec, x_data, y_data, num_features, int_feats, 
                           color_map, box_cat, box_cols){
  rf_imp <- data.frame(Overall = varImp(mod)$importance$Overall)
  
  good_gene_tab <- feat_vec

  for(b in 1:ncol(x_data)){
    if(grepl(';',colnames(x_data)[b])){
      lookerupper <- paste0("^",colnames(x_data)[b], "$")
      if(grepl("\\[", lookerupper)){
        lookerupper1 <- strsplit(lookerupper, "\\[")
        lookerupper2 <- strsplit(lookerupper, "\\]")
        lookerupper <- paste0(lookerupper1[[1]][1], "\\[.*\\]",
                              lookerupper2[[1]][2])
      }
      lookupnum <- grep(lookerupper, good_gene_tab$Feature, fixed=F)
      colnames(x_data)[b] <- paste0("OTU_", lookupnum)
    }
  }
  feat_heat <- data.frame(Features = colnames(x_data))
  feat_heat <- cbind(feat_heat, rf_imp)
  feat_heat$DataType = c(rep("Gene", ncol(x_data_genes)),
                         rep("OTU", ncol(x_data_otus)),
                         rep("ImmunePop", ncol(x_data_facs)),
                         rep("Cytokine", ncol(x_data_cyt)))
  
  feat_heat_filt <- feat_heat[order(feat_heat$Overall, decreasing=T),][1:num_features,]
  feat_heat_filt$Features <- factor(feat_heat_filt$Features, levels = rev(feat_heat_filt$Features))
  
  #pdf("All_feature_imp_clean.pdf", height = 4, width = 6)
  g_imp = ggplot(feat_heat_filt, aes(Features, Overall, fill=DataType)) +
    geom_col() + coord_flip() + 
    theme_bw() +
    ylab("MeanDecreaseGini") + ggtitle(paste0("Top ", num_features, " Model Features")) +
    scale_fill_manual(values=color_map)+
    theme(axis.text.y = element_text(size=12, color="black"),#, angle = 60, hjust=1),
          axis.text.x = element_blank(),
          axis.title = element_text(size=13, color="black"),
          legend.text = element_text(size=12, color="black"),
          legend.title = element_text(size=13, color="black"),
          plot.title = element_text(size=14, face='bold'))
  #dev.off()
  
  
  #### imp boxes
  
  x_small <- x_data[,as.character(feat_heat_filt$Features)]
  
  box_df <- data.frame(x_small[,int_feats], cat=y_data)
  box_df <- melt(box_df)
  g=ggplot(box_df, aes(cat, log2(value+1), color=cat, group=cat)) +
    geom_boxplot(alpha=0.5, outlier.shape=NA) +
    geom_jitter(width = 0.2) +
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
  box_exp = g+facet_wrap(~variable, nrow = 1, scales='free') +
    theme(strip.text = element_text(size=15))
  
  ####network
  
  
 
  top_feats <- as.character(feat_heat[order(feat_heat$Overall, decreasing=T),][1:num_features,]$Features)
  top_feats_meta <- feat_heat[order(feat_heat$Overall, decreasing=T),][1:num_features,]
  rownames(top_feats_meta) <- top_feats_meta$Features
  
  x_data_top <- x_data[,top_feats]
  mat=cor(x_data_top)
  mat[mat<0.6 & mat > -0.6]=0
  good_cols <- which(colSums(abs(mat))>1)
  
  mat=mat[good_cols,good_cols]
  
  ## edge thickness
  connecs = c()
  for(i in 1:ncol(mat)){
    huh = length(which(mat[,i]!=0))
    connecs <- c(connecs, huh)
  }
  
  p_plot <- data.frame(names=rownames(mat))
  p_plot$edges = connecs
  p_plot$edge_adj <- (round(p_plot$edges/1,0)+1)*3
  
  
  #
  
  network=graph_from_adjacency_matrix(as.matrix(mat), weighted=T, mode="undirected", diag=F)
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
  
  coul <- c(Cytokine="darkseagreen2", Gene="lightcoral", ImmunePop="dodgerblue2", OTU="mediumorchid3")
  my_color=coul[top_feats_meta$DataType[good_cols]]
  # plot
  
  curr_coord = icoord
  new_coord = data.frame(curr_coord, top_feats_meta[good_cols,])
  colnames(new_coord) <- c("X", "Y", 'Feature',"Importance", "DataType")
  
  #library(ggplot2)
  #ggplot(new_coord, aes(X,Y,color=DataType)) + 
  #  geom_point()
  
  gs = as.character(unique(new_coord$DataType))
  gs=c("Gene", "ImmunePop", "Cytokine", "OTU")
  better_coords <- data.frame(X=NA,Y=NA,Feature=NA,Importance=NA,DataType=NA,X2=NA,Y2=NA)
  low_bounds_y=c(250,800,250,0)
  high_bounds_y=c(700,950,700,200)
  
  low_bounds_x=c(0,100,700,300)
  high_bounds_x=c(300,700,1000,700)
  
  for(i in 1:length(gs)){
    try({
    curr <- subset(new_coord, DataType == gs[i])
    
    x2 = sample(low_bounds_x[i]:high_bounds_x[i], nrow(curr))
    y2 = sample(low_bounds_y[i]:high_bounds_y[i], nrow(curr), replace = T)
    
    curr$X2=x2
    curr$Y2=y2
    
    better_coords = rbind(better_coords, curr)
    })
  }
  better_coords=better_coords[-1,]
  
  #ggplot(better_coords, aes(X2,Y2,color=DataType)) + 
  #  geom_point()
  better_coords=better_coords[rownames(new_coord),]
  
  return(list(g_imp, box_exp, network, better_coords, my_color, p_plot$edge_adj))
  
  #pdf("int_network5.26.pdf", height = 50, width = 50)
  #par(bg="white", mar=c(0,0,0,0))
  #set.seed(4)
  #net_net=as.grob(
  #  function()
  #    plot(network, layout=as.matrix(better_coords[,6:7]),
  #         vertex.size=10,
  #         vertex.color=my_color, 
  #         vertex.label.cex=1.2,
  #         vertex.label.family="Helvetica",
  #         vertex.label.color="black",
  #         vertex.frame.color="transparent"
  #    )
  #  )
}

meta <- readRDS("ml_inputs/mouse_metadata.RData")
y_data <- meta$Environment

x_data_genes <- readRDS("ml_inputs/genes.RData")
# By Variation
x_data_genes <- x_data_genes[,order(colVars(x_data_genes),decreasing=T)][,1:200]

##otus
x_data_otus <- readRDS("ml_inputs/otus.RData")
##facs
x_data_facs <- readRDS("ml_inputs/facs.RData")
##cytokines
x_data_cyt <- readRDS("ml_inputs/cytokines.RData")

#### all
x_data_all <- cbind(x_data_genes, x_data_otus, x_data_facs, x_data_cyt)


mod <- readRDS("predict_env/All_rf_model.RData")
feat_vec <- read.table("predict_env/All_feature_vector.txt", T, '\t')
num_features = 40
color_map <- c(Cytokine="darkseagreen2", Gene="lightcoral", ImmunePop="dodgerblue2", OTU="mediumorchid3")

int_feats <- c("MLN_CD8_T_cells", "MLN_CD4_T_cells", "Nfkbia", "Egr1", "Atf4")

predict_env <- mod_importance(mod,feat_vec,num_features = 40, int_feats=int_feats, 
                              x_data=x_data_all, 
                              y_data=y_data, color_map=color_map, box_cat = "Environment", 
                              box_cols = c("lab"="mediumorchid3", "wild"="red3"))

tkplot(predict_env[[3]], layout=as.matrix(predict_env[[4]][,6:7]), 
       vertex.size=predict_env[[6]],
       vertex.color=predict_env[[5]], 
       vertex.label.cex=1.3,
       vertex.label.family="Helvetica",
       vertex.label.color="black"#,
       #vertex.frame.color="transparent"
)

tk_coords = tk_coords(2)

predict_env[[3]] = set_edge_attr(predict_env[[3]], "label", value="")
net_net <- as.grob(function()
  plot(predict_env[[3]], #layout=as.matrix(better_coords[,5:6]),
     layout=tk_coords,
     vertex.size=predict_env[[6]],
     vertex.color=predict_env[[5]], 
     vertex.label.cex=1.1,
     vertex.label.family="Helvetica",
     vertex.label.color="black",
     vertex.frame.color="transparent"
)
)

pdf("env_model.pdf", height = 8, width = 17)
g1=arrangeGrob(net_net, predict_env[[2]], nrow=2)
grid.arrange(predict_env[[1]], g1, nrow=1, widths = c(1,1.8))
dev.off()

pdf("env_model_7.19.pdf", height = 8, width = 17)
grid.arrange(predict_env[[1]], net_net, nrow=1, widths=c(1,1.5))
#grid.arrange(predict_env[[1]], g1, nrow=1, widths = c(1,1.8))
dev.off()

#### stats

int_things <- c("MLN_CD8_T_cells", "MLN_CD4_T_cells", "Nfkbia", "Egr1", "Atf4")
stats = data.frame(p=NA)
for(i in 1:length(int_things)){
  x_sub <- x_data_all[,int_things[i]]
  pp=wilcox.test(x_sub ~ meta$Environment)$p.value
  stats <- rbind(stats, data.frame(p=pp))
}
stats = stats[-1,]


#

#
###
meta <- readRDS("ml_inputs/mouse_metadata.RData")
y_data <- meta$Genotype

mod <- readRDS("predict_genotype//All_rf_model.RData")
feat_vec <- read.table("predict_genotype//All_feature_vector.txt", T, '\t')
num_features = 40
color_map <- c(Cytokine="darkseagreen2", Gene="lightcoral", ImmunePop="dodgerblue2", OTU="mediumorchid3")

int_feats <- c("MLN_CD8_T_cells_CD44hi_CD62Lhi", "MLN_Ly6C_Monocyte", 
               "Erdr1", "IL.10_PseudomonasA", "IL.5_PBS")

predict_geno <- mod_importance(mod,feat_vec,num_features = 40, int_feats=int_feats, 
                              x_data=x_data_all, 
                              y_data=y_data, color_map=color_map, box_cat = "Genotype", 
                              box_cols = c("AtgE"="dodgerblue2", "AtgH"="navyblue",
                                           "B6"="darkorange1", "NOD2"="mediumspringgreen"))

tkplot(predict_geno[[3]], layout=as.matrix(predict_geno[[4]][,6:7]), 
       vertex.size=predict_geno[[6]],
       vertex.color=predict_geno[[5]], 
       vertex.label.cex=1.3,
       vertex.label.family="Helvetica",
       vertex.label.color="black"#,
       #vertex.frame.color="transparent"
)

tk_coords = tk_coords(5)

predict_geno[[3]] = set_edge_attr(predict_geno[[3]], "label", value="")
net_net <- as.grob(function()
  plot(predict_geno[[3]], #layout=as.matrix(better_coords[,5:6]),
       layout=tk_coords,
       vertex.size=predict_geno[[6]],
       vertex.color=predict_geno[[5]], 
       vertex.label.cex=1.1,
       vertex.label.family="Helvetica",
       vertex.label.color="black",
       vertex.frame.color="transparent"
  )
)

pdf("geno_model.pdf", height = 8, width = 18)
g1=arrangeGrob(net_net, predict_geno[[2]], nrow=2)
grid.arrange(predict_geno[[1]], g1, nrow=1, widths = c(1,1.8))
dev.off()

pdf("geno_model_7.19.pdf", height = 8, width = 17)
grid.arrange(predict_geno[[1]], net_net, nrow=1, widths=c(1,1.5))
dev.off()

int_things <- c("MLN_CD8_T_cells_CD44hi_CD62Lhi", "MLN_Ly6C_Monocyte", "Erdr1", "IL.10_PseudomonasA", "IL.5_PBS")
stats = data.frame(p=NA)
for(i in 1:length(int_things)){
  x_sub <- x_data_all[,int_things[i]]
  pp=aov(x_sub ~ meta$Genotype)
  p=summary(pp)[[1]][["Pr(>F)"]][1]
  
  if(p<0.05){
  thing = data.frame(TukeyHSD(pp)$`meta$Genotype`)
  print(i)
  print(thing)
  }
}
stats = stats[-1,]
