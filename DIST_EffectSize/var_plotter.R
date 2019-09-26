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


type_cols <- c(Diarrhea="red2", Environment="dodgerblue2", FACS = "forestgreen",
               Wedge_cage="navyblue",Flow.date="grey50",
               Gender="hotpink",Genotype="darkorange1", Pregnant="purple", WeightGain="forestgreen")


#function for plotting PCoA with Effect Sizes
var_plotter <- function(input, effectors, col_var, col_var_name="Environment",
                        shape_var=NA, shape_var_name=NA, Feature, dist_panel=F){
  
  D <- dist(input, method = "euclidean")
  res <- pcoa(D)
  plotter <- data.frame(res$vectors, col_var=col_var, 
                        shape_var = shape_var)
  
  #good_plotter <- subset(plotter, Axis.1 > quantile(plotter$Axis.1, 0.05) &
  #                         Axis.1 < quantile(plotter$Axis.1, 0.95))
  #good_plotter <- subset(good_plotter, Axis.2 > quantile(plotter$Axis.2, 0.05) &
  #                         Axis.2 < quantile(plotter$Axis.2, 0.95))

  
  g0=ggplot(plotter, aes(Axis.1, Axis.2, color=col_var, shape=shape_var)) + 
    geom_point(size=5) +
    guides(color=guide_legend(title=col_var_name), 
           shape=guide_legend(title=shape_var_name)) +
    ylab(paste0("PCo2 ", round(res$values$Broken_stick[2]*100,2), "% expl. variation")) +
    xlab(paste0("PCo1 ", round(res$values$Broken_stick[1]*100,2), "% expl. variation")) +
    scale_color_manual(values = c("mediumpurple1", "red3","dodgerblue2")) +
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
    scale_color_manual(values = c("mediumpurple1", "red3","dodgerblue2")) +
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
    scale_color_manual(values = c("red3", "mediumorchid3")) + coord_equal() + 
    geom_segment(data=datapc2, aes(x=0, y=0, xend=Axis.1, yend=Axis.2), 
                 arrow=arrow(length=unit(0.2,"cm")), alpha=0.5) +
    geom_label_repel(data=datapc2, aes(x=Axis.1, y=Axis.2, label=varnames), 
                     size = 4, force=4, segment.alpha=0.5) +
    theme_void()
  #g_seg
  
  
  
  pc1_hist <- ggplot(plotter, aes(Axis.1, color = col_var, fill= col_var)) +
    geom_density(alpha=0.6) + 
    scale_y_reverse() +
    scale_color_manual(values = c("mediumpurple1", "red3","dodgerblue2")) +
    scale_fill_manual(values = c("mediumpurple1", "red3","dodgerblue2")) +
    theme_void()
  
  pc2_hist <- ggplot(plotter, aes(Axis.2, color = col_var, fill= col_var)) +
    geom_density(alpha=0.6) + 
    coord_flip() +
    scale_y_reverse() +
    scale_color_manual(values = c("mediumpurple1", "red3", "dodgerblue2")) +
    scale_fill_manual(values = c("mediumpurple1", "red3", "dodgerblue2")) +
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

  
#Effect size calculation  
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
    
    return(arrangeGrob(pcoa_plot, g_seg_plot,dist_plot, gg_combo, ncol=4, widths=c(2,2,2,1.7)))
    
  } else {
    return(arrangeGrob(pcoa_plot, g_seg_plot, gg_combo, ncol=3, widths=c(2,2,1.5)))
  }
  
}
