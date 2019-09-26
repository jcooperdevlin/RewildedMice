#### MLN Stim Fix 5.9

setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/data/MLN_stimulations/mln_fix")

library(xlsx)
library(ggplot2)
library(drc)
library(drLumi)
library(nplr)
library(ELISAtools)
library(circlize)
library(cowplot)
library(gridExtra)
library(ggsci)
library(reshape)
library(ggrepel)

# read in data (wrtie a function that does this for each plate?)

cba_reader <- function(plate_files=NA, analyte_indices=NA, plate_names=NA){
  sample_data <- data.frame(plate=NA, well=NA, analyte=NA, sample=NA, mfi=NA)
  standard_data <- data.frame(sample=NA, analyte=NA, plate=NA, ec=NA, mfi=NA)
  
  for(j in 1:length(plate_files)){
    print(paste0("Working on ", plate_names[j]))
    for (i in analyte_indices){
      name <- colnames(read.xlsx(plate_files[j], sheetIndex = i, startRow = 0, endRow = 22))[1]
      print(paste0("Reading in ", name, "on ", plate_names[j]))
      st <- read.xlsx(plate_files[j], sheetIndex = i, startRow = 5, endRow = 22)
      st$plate <- rep(plate_names[j], nrow(st))
      st$analyte <- rep(name, nrow(st))
      exp <- unique(st$Expected.pg.ml.)
      st$ec <- rep(exp[c(1,3:length(exp))], each = 2)
      st$well = NA
      st$sample <- rep(paste0("Standard", 1:(nrow(st)/2)),each=2)
      st$mfi <- st$MFI
      
      
      stand_add <- st[,c("sample", "analyte", "plate", "ec", "mfi")]
      standard_data <- rbind(standard_data, stand_add)
      
      
      ## samples
      name <- colnames(read.xlsx(plate_files[j], sheetIndex = i, startRow = 0, endRow = 22))[1]
      samp <- read.xlsx(plate_files[j], sheetIndex = i, startRow = 24)
      samp$plate <- rep(plate_names[j], nrow(samp))
      samp$analyte <- rep(name, nrow(samp))
      samp$well = NA
      samp$sample <- paste0(samp$Mice, ";", samp$Stimulation)
      samp$sample <- gsub(" ", "", samp$sample)
      samp$mfi <- samp$MFI
      samp <- samp[!is.na(samp$Name),]
      
      
      samp_add <- samp[,c("plate", "well", "analyte", "sample", "mfi")]
      sample_data <- rbind(sample_data, samp_add)
    }
  }
  standard_data <- standard_data[-1,]
  sample_data <- sample_data[-1,]
  
  return(list(standard_data, sample_data))
}

#
#
#
#
# read in serum data

plate_files <- list.files("mln_input", pattern="Plate", full.names = T)
plate_files <- plate_files[grepl("~", plate_files)==F]
analyte_indices <- c(4:16)
plate_names <- c("Plate-11", "Plate-12", "Plate-13", "Plate-14", "Plate-5", 
                 "Plate-6", "Plate-7", "Plate-8","Plate-9")

serum_data = cba_reader(plate_files = plate_files, 
                        analyte_indices = analyte_indices, 
                        plate_names=plate_names)



#
#
#

fitter_combine_plates <- function(analyte_name="Analyte", input_mfi=NA, input_ec=NA){
  plist=list()
  colors_clusters = c(pal_d3("category10")(10), pal_d3("category20b")(20), pal_igv("default")(51))
  og_data = ggplot(input_ec, aes(ec+1, mfi+1, color=plate)) +
    geom_point() + scale_x_log10() + scale_y_log10() +
    xlab("Log10 Concentration") + ylab("Log10 MFI") +
    scale_color_manual(values=colors_clusters) +
    theme_bw()
  
  #
  # fit the data
  #stand_df <- subset(input_ec, plate == plates[i] & !is.na(mfi))
  #stand_df <- aggregate(mfi ~ sample + analyte + ec, mean, data=stand_df)
  
  stand_df <- subset(input_ec, !is.na(mfi))
  
  np_fit = nplr(log10(stand_df$ec+1), log10(stand_df$mfi+1), useLog = F)
  summary(np_fit)
  
  x <- getXcurve(np_fit)
  y <-getYcurve(np_fit)
  
  newdata <- data.frame(x=x,y=y)
  
  #calculate residual by hand
  res <- data.frame(
    xs = np_fit@x,
    y_act = np_fit@y,
    y_fit = np_fit@yFit
  )
  res$res <- abs(res$y_act - res$y_fit)
  res$color = "black"
  res$color[res$res>0.1] <- "red"
  res$color = factor(res$color, levels = c("black", "red"))
  
  g_res <- ggplot(res, aes(y_act, res, color = color)) + 
    geom_point(show.legend = F) + coord_flip() +
    geom_hline(yintercept = 0.1) +
    ggtitle("Original Residuals") +
    ylab("residuals") + xlab("log10 MFI") +
    scale_color_manual(values=c("black", "red")) +
    theme_bw()
  
  stand_df <- stand_df[order(stand_df$ec, stand_df$plate),]
  stand_df$point_color = res$color
  stand_df$point_color = factor(stand_df$point_color, levels = c("black", "red"))
  g_fit=ggplot(stand_df, aes(log10(ec+1), log10(mfi+1), color=point_color)) +
    geom_point(show.legend = F) + 
    scale_color_manual(values=c("black", "red")) +
    ggtitle("Original Fit") +
    geom_line(aes(x=x,y=y), newdata, inherit.aes = F) +
    xlab("Log10 Concentration") + ylab("Log10 MFI") +
    theme_bw()
  
  fit_plot1=arrangeGrob(g_fit, g_res, nrow=1)
  #
  # based on this plot ^^^
  # if points fall off standrd curve remove them and refit
  
  if(any(stand_df$point_color=="red")) {
    stand_df <- subset(stand_df, point_color != "red")
    #stand_df <- aggregate(mfi ~ sample + analyte + ec, mean, data=stand_df)
    
    np_fit = nplr(log10(stand_df$ec+1), log10(stand_df$mfi+1), useLog = F)
    summary(np_fit)
    
    x <- getXcurve(np_fit)
    y <-getYcurve(np_fit)
    
    newdata <- data.frame(x=x,y=y)
    
    #calculate residual by hand
    res <- data.frame(
      xs = np_fit@x,
      y_act = np_fit@y,
      y_fit = np_fit@yFit
    )
    res$res <- abs(res$y_act - res$y_fit)
    res$color = "black"
    res$color[res$res>0.1] <- "red"
    res$color = factor(res$color, levels = c("black", "red"))
    
    g_res <- ggplot(res, aes(y_act, res, color = color)) + 
      geom_point(show.legend = F) + coord_flip() +
      geom_hline(yintercept = 0.1) +
      ggtitle("Corrected Residuals") +
      ylab("residuals") + xlab("log10 MFI") +
      scale_color_manual(values=c("black", "red")) +
      theme_bw()
    
    stand_df <- stand_df[order(stand_df$ec, stand_df$plate),]
    stand_df$point_color = res$color
    stand_df$point_color = factor(stand_df$point_color, levels = c("black", "red"))
    g_fit=ggplot(stand_df, aes(log10(ec+1), log10(mfi+1), color=point_color)) +
      geom_point(show.legend = F) + 
      ggtitle("Corrected Fit") +
      scale_color_manual(values=c("black", "red")) +
      geom_line(aes(x=x,y=y), newdata, inherit.aes = F) +
      xlab("Log10 Concentration") + ylab("Log10 MFI") +
      theme_bw()
    
    fit_plot2=arrangeGrob(g_fit, g_res, nrow=1)
    
  } else {fit_plot2=ggplot()+theme_void()}
  
  
  new_points = getEstimates(np_fit, log10(stand_df$mfi+1))
  
  new_points$act_c <- stand_df$ec
  new_points$est_c <- 10^new_points$x
  new_points$recovery <- (new_points$act_c/new_points$est_c)*100
  
  myPalette <- colorRampPalette("dodgerblue2")
  r_plot=ggplot(stand_df, aes(log10(ec+1), log10(mfi+1))) +
    geom_point() + 
    ggtitle("Standard Recovery") +
    geom_line(aes(x=x,y=y), newdata) +
    geom_point(aes(x=x,y=y, col=recovery), new_points, show.legend = F) +
    scale_colour_gradientn(colours = myPalette(n=10), limits=c(80, 120)) +
    xlab("Log10 Concentration") + ylab("Log10 MFI") +
    theme_bw()
  
  #
  #
  #fit unknowns and plot according to fit
  
  sample_df <- subset(input_mfi, !is.na(mfi))
  sample_points = getEstimates(np_fit, log10(sample_df$mfi+1))
  
  est_c <- 10^sample_points$x
  sample_points$est_c <- est_c
  
  s_plot=ggplot(sample_points, aes(x, y)) +
    geom_point() + 
    ggtitle("Samples projected on Fit") +
    geom_line(aes(x=x,y=y), newdata, alpha=0.5) +
    xlab("Log10 Concentration") + ylab("Log10 MFI") +
    theme_bw()
  
  plist[[1]]=arrangeGrob(og_data,
                         arrangeGrob(fit_plot1, fit_plot2, 
                                     arrangeGrob(r_plot, s_plot,nrow=1), 
                                     nrow = 3),
                         ncol=2)
  return(list(plist, est_c))
}
#
#
### run through the analytes


#plate combine
analytes <- unique(serum_data[[2]]$analyte)
serum_fix <- subset(serum_data[[2]], analyte == analytes[1])[,1:4]
for(k in 1:length(analytes)){
  analyte_name = analytes[k]
  input_mfi <- subset(serum_data[[2]], analyte == analyte_name & !is.na(mfi))
  input_ec <- subset(serum_data[[1]], analyte == analyte_name)
  
  full_plots <- fitter_combine_plates(analyte_name, input_mfi, input_ec)
  namer=paste0("plots_combine/", analyte_name, "_fit-correction.pdf")
  pdf(namer, height = 7, width = 12)
  pp=plot_grid(plotlist = full_plots[[1]], ncol=1)
  print(pp)
  dev.off()
  
  input_mfi$est_c <- full_plots[[2]]
  rownames(input_mfi) <- as.character(input_mfi$sample)
  input_mfi <- input_mfi[as.character(serum_fix$sample),]
  #serum_fix2 <- merge(serum_fix, input_mfi[,c("sample", "est_c")], by="sample", all.x=T)
  serum_fix <- cbind(serum_fix, input_mfi$est_c)
  colnames(serum_fix)[ncol(serum_fix)] <- analytes[k]
}


#
#
#
#
#
#
#



#melt and cast right variables
stim_melt <- melt(serum_fix, measure.vars = colnames(serum_fix)[5:17], id.vars = colnames(serum_fix)[c(1,4)])
stim_melt$value[is.na(stim_melt$value)] <- 0
stim_melt$challenge <- gsub(".*;", "", stim_melt$sample)
stim_melt$mouse_id <- gsub(";.*", "", stim_melt$sample)
stim_melt$variable <- gsub("Analyte\\.", "", stim_melt$variable)
stim_cast <- cast(stim_melt, mouse_id + plate ~ variable + challenge, mean)

stim_cast$mouse_id <- gsub(" ", "", stim_cast$mouse_id)

#
#
rownames(stim_cast) <- stim_cast$mouse_id



meta <- read.table("../../../../int/data/metadata/mice_metadata.11.19_mouse_id.txt", T, '\t')
rownames(meta) <- as.character(meta$mouse_id)
cyt_keep <- intersect(stim_cast$mouse_id, meta$mouse_id)


mln_cyt_keep <- stim_cast[cyt_keep,]
mln_cyt_keep[is.na(mln_cyt_keep)] <- 1
mln_cyt_keep[mln_cyt_keep<=1] <- 0
meta_keep <- meta[cyt_keep,]

good_cols <- c(3:106)
mln_cyt_keep = mln_cyt_keep[,c(1,good_cols)]

write.table(mln_cyt_keep,"MLN_names_final.txt", sep='\t', quote = F, row.names=F)

#

#

pc <- prcomp(log2(mln_cyt_keep[,-1]+1), scale=T)
pc <- prcomp(log2(mln_cyt_keep[,3:106]+1), scale=T)
huh<-pc$rotation

plotter <- data.frame(PC1=pc$x[,1], PC2=pc$x[,2], mln_cyt_keep, meta_keep)

ggplot(plotter, aes(PC1, PC2, color = Environment)) +
  geom_point()

ggplot(plotter, aes(PC1, PC2, color = plate)) +
  geom_point()

ggplot(plotter, aes(PC1, PC2, color = factor(Wedge_cage))) +
  geom_point()

ggplot(plotter, aes(PC1, PC2, color = Genotype)) +
  geom_point()

ggplot(plotter, aes(PC1, PC2, color = Gender)) +
  geom_point()

table(meta_keep$Gender, meta_keep$Environment)
table(meta_keep$Genotype, mln_cyt_keep$plate)


library(ComplexHeatmap)

pdf("mln_plate_gender_bias.pdf", height = 8, width = 15)
va <- rowAnnotation(df = data.frame(Plate=mln_cyt_keep$plate, Gender=meta_keep$Gender))
Heatmap(log2(mln_cyt_keep[,3:106]+1), 
        col = colorRamp2(c(1, 7), c("white", "red")),
        cluster_columns = T, cluster_rows = T,
        show_row_names = T, show_column_names = T,
        row_names_gp = gpar(fontsize = 6),
        heatmap_legend_param=list(title = "Log2 Value")) + va
dev.off()

pdf("mln_plate_env.pdf", height = 8, width = 15)
va <- rowAnnotation(df = data.frame(Environment=plotter$Environment, Genotype=plotter$Genotype))
Heatmap(log2(mln_cyt_keep[,3:106]+1), 
        col = colorRamp2(c(1, 7), c("white", "red")),
        cluster_columns = T, cluster_rows = T,
        show_row_names = T, show_column_names = T,
        row_names_gp = gpar(fontsize = 6),
        heatmap_legend_param=list(title = "Log2 Value")) + va
dev.off()


#
#
#
#
#

#
#
#####
#save time read it in
mln_cyt_keep = read.table("MLN_names_final.txt", T, '\t')

meta <- read.table("../../../../int/data/metadata/mice_metadata.11.19_mouse_id.txt", T, '\t')
rownames(meta) <- as.character(meta$mouse_id)
meta_keep <- meta[as.character(mln_cyt_keep$mouse_id),]

pc <- prcomp(log2(mln_cyt_keep[,-1]+1), scale=T)
huh<-pc$rotation

plotter <- data.frame(PC1=pc$x[,1], PC2=pc$x[,2], mln_cyt_keep, meta_keep)

##
stim_melt <- melt(plotter, measure.vars = colnames(plotter)[4:107], id.vars = colnames(plotter)[c(3,109,110,113)])
stim_melt$challenge <- gsub(".*_", "", stim_melt$variable)
stim_melt$cytokine <- gsub("_.*", "", stim_melt$variable)
stim_melt$cytokine <- gsub(".*MCP", "MCP", stim_melt$cytokine)
stim_melt$cytokine <- gsub(".*MIP", "MIP", stim_melt$cytokine)
stim_melt$cytokine <- gsub(".*IL", "IL", stim_melt$cytokine)
stim_melt$cytokine <- gsub(".*KC", "KC", stim_melt$cytokine)
stim_melt$cytokine <- gsub(".*IFN", "IFN", stim_melt$cytokine)
stim_melt$cytokine <- gsub(".*TNF", "TNF", stim_melt$cytokine)

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
stim_melt$cytokine <- gsub(".*MCP", "MCP", stim_melt$cytokine)
stim_melt$cytokine <- gsub(".*MIP", "MIP", stim_melt$cytokine)
stim_melt$cytokine <- gsub(".*IL", "IL", stim_melt$cytokine)
stim_melt$cytokine <- gsub(".*KC", "KC", stim_melt$cytokine)
stim_melt$cytokine <- gsub(".*IFN", "IFN", stim_melt$cytokine)
stim_melt$cytokine <- gsub(".*TNF", "TNF", stim_melt$cytokine)

challenges <- c("BacillusS","BacteroidesV","CandidaA",
                "CD3.CD28","ClostridiumP","PseudomonasA", 
                "StaphA", "PBS")
glist = list()
for (j in 1:length(challenges)){
  full_df2 <- subset(stim_melt, challenge == challenges[j])
  
  nums2 <- log2(full_df2$value+1)
  
  stuff <- full_df2$Environment
  cyto <- full_df2$cytokine
  
  plotter <- data.frame(stuff, cyto, nums2)
  
  plotter_slim <- aggregate(plotter$nums2, list(plotter$stuff, plotter$cyto),mean)
  colnames(plotter_slim) <- c("Environment", "cyto", 'nums')
  
  plotter_slim$Environment <- factor(plotter_slim$Environment, levels = c("lab", "wild"))
  plotter_slim$cyto <- gsub("\\.\\.", "_", plotter_slim$cyto)
  plotter_slim$cyto <- gsub("\\.", "_", plotter_slim$cyto)
  plotter_slim$cyto <- substr(plotter_slim$cyto,1,nchar(plotter_slim$cyto)-1)
  plotter_slim=plotter_slim[order(plotter_slim$nums, decreasing = T),]
  plotter_slim$cyto <- factor(plotter_slim$cyto, levels = rev(unique(plotter_slim$cyto)))
  plotter_slim=plotter_slim[order(plotter_slim$cyto, decreasing = T),]
  plotter_slim$labelery <- rep(max(plotter_slim$nums), nrow(plotter_slim))
  plotter_slim$labelerx <- ""
  for(i in 1:nrow(plotter_slim)){
    if(plotter_slim$Environment[i]=="wild"){
      plotter_slim$labelerx[i]<-as.character(plotter_slim$cyto[i])
    }
  }
  
  #filenamer <- paste0("radars/radar_plot_", challenges[j], ".pdf")
  #pdf(filenamer)
  glist[[j]] <- ggplot(plotter_slim, aes(x=cyto, y=nums)) +
    geom_polygon(aes(group = Environment, color = Environment,fill = Environment),
                 alpha=0.4, size = 1, show.legend = TRUE) +
    coord_polar(clip='on')+
    geom_text(data=plotter_slim, aes(x=cyto, y=labelery, label=labelerx), size = 5) + 
    scale_color_manual(values = c("mediumpurple1", "red3"))+
    scale_fill_manual(values = c("mediumpurple1", "red3"))+
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


pdf("all_radars.pdf", height = 12, width = 32)
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
  plotter_slim$cyto <- gsub("\\.\\.", "_", plotter_slim$cyto)
  plotter_slim$cyto <- gsub("\\.", "_", plotter_slim$cyto)
  plotter_slim$cyto <- substr(plotter_slim$cyto,1,nchar(plotter_slim$cyto)-1)
  plotter_slim=plotter_slim[order(plotter_slim$nums, decreasing = T),]
  plotter_slim$cyto <- factor(plotter_slim$cyto, levels = rev(unique(plotter_slim$cyto)))
  plotter_slim=plotter_slim[order(plotter_slim$cyto, decreasing = T),]
  plotter_slim$labelery <- rep(max(plotter_slim$nums), nrow(plotter_slim))
  
  
  #filenamer <- paste0("radars/radar_plot_", challenges[j], ".pdf")
  #pdf(filenamer)
  glist[[j]] <- ggplot(plotter_slim, aes(x=cyto, y=nums)) +
    geom_polygon(aes(group = Genotype, color = Genotype,fill = Genotype),
                 alpha=0.4, size = 1, show.legend = TRUE) +
    coord_polar(clip='on')+
    geom_text(data=plotter_slim, aes(x=cyto, y=labelery, label=cyto), size = 4) + 
    #scale_color_manual(values = c("mediumpurple1", "red3"))+
    #scale_fill_manual(values = c("mediumpurple1", "red3"))+
    xlab("Cytokines") + ylab("Log pg/mL") +
    ggtitle(challenges[j]) +
    theme(
      #plot.margin=margin(10, 10, 5, 5),
      panel.grid.major.x = element_line(colour = "grey90"),
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
