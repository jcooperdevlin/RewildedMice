#### Serum Fix 4.29

###
setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/data/Serum_stimulations/serum_fix")

library(xlsx)
library(ggplot2)
library(drc)
library(drLumi)
library(nplr)
library(ELISAtools)
library(circlize)
library(cowplot)
library(gridExtra)

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
      samp$sample <- gsub(" ", "", samp$Name)
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

plate_files <- list.files(".", pattern="SerumCBA", full.names = T)
plate_files <- plate_files[grepl("~", plate_files)==F]
analyte_indices <- c(4:16)
plate_names <- c("Plate-1", "Plate-2")

serum_data = cba_reader(plate_files = plate_files, 
                        analyte_indices = analyte_indices, 
                        plate_names=plate_names)


#
#
#
#
#
#
#
#write a function that 
#   1. fits standard curve, 
#   2. evaluates fit
#   3. removes points if necessary
#   4. Refits

#
#
#


fitter <- function(analyte_name="Analyte", input_mfi=NA, input_ec=NA){
  
  og_data = ggplot(input_ec, aes(ec+1, mfi+1, color=plate)) +
    geom_point() + scale_x_log10() + scale_y_log10() +
    xlab("Log10 Concentration") + ylab("Log10 MFI") +
    scale_color_manual(values=c("purple", "dodgerblue2")) +
    theme_bw()
  
  #
  # fit the data
  plates <- unique(input_ec$plate)
  plist=list()
  for(i in 1:length(plates)){
    stand_df <- subset(input_ec, plate == plates[i] & !is.na(mfi))
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
    res$color[res$res>0.05] <- "red"
    res$color = factor(res$color, levels = c("black", "red"))
    
    g_res <- ggplot(res, aes(y_act, res, color = color)) + 
      geom_point(show.legend = F) + coord_flip() +
      geom_hline(yintercept = 0.05) +
      ggtitle("Original Residuals") +
      ylab("residuals") + xlab("log10 MFI") +
      scale_color_manual(values=c("black", "red")) +
      theme_bw()
    
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
      res$color[res$res>0.05] <- "red"
      res$color = factor(res$color, levels = c("black", "red"))
      
      g_res <- ggplot(res, aes(y_act, res, color = color)) + 
        geom_point(show.legend = F) + coord_flip() +
        geom_hline(yintercept = 0.05) +
        ggtitle("Corrected Residuals") +
        ylab("residuals") + xlab("log10 MFI") +
        scale_color_manual(values=c("black", "red")) +
        theme_bw()
      
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
    
    sample_df <- subset(input_mfi, plate == plates[i] & !is.na(mfi))
    sample_points = getEstimates(np_fit, log10(sample_df$mfi+1))
    
    sample_points$est_c <- 10^sample_points$x

    s_plot=ggplot(sample_points, aes(x, y)) +
      geom_point() + 
      ggtitle("Samples projected on Fit") +
      geom_line(aes(x=x,y=y), newdata, alpha=0.5) +
      xlab("Log10 Concentration") + ylab("Log10 MFI") +
      theme_bw()
    
    plist[[i]]=arrangeGrob(fit_plot1, fit_plot2, 
                           arrangeGrob(r_plot, s_plot,nrow=1), 
                           nrow = 3, top = plates[i])
  }
  return(plist)
}


#
#

#### 
#
#
#
fitter_combine_plates <- function(analyte_name="Analyte", input_mfi=NA, input_ec=NA){
  plist=list()
  og_data = ggplot(input_ec, aes(ec+1, mfi+1, color=plate)) +
    geom_point() + scale_x_log10() + scale_y_log10() +
    xlab("Log10 Concentration") + ylab("Log10 MFI") +
    scale_color_manual(values=c("purple", "dodgerblue2")) +
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

#analytes <- unique(serum_data[[2]]$analyte)
#for(k in 1:length(analytes)){
#  analyte_name = analytes[k]
#  input_mfi <- subset(serum_data[[2]], analyte == analyte_name)
#  input_ec <- subset(serum_data[[1]], analyte == analyte_name)
#  
#  full_plots <- fitter(analyte_name, input_mfi, input_ec)
#  namer=paste0("plots/", analyte_name, "_fit-correction.pdf")
#  pdf(namer, height = 7, width = 12)
#  pp=plot_grid(plotlist = full_plots, ncol=2)
#  print(pp)
#  dev.off()
#}

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
  serum_fix2 <- merge(serum_fix, input_mfi[,c("sample", "est_c")], by="sample", all.x=T)
  serum_fix <- cbind(serum_fix, serum_fix2[,ncol(serum_fix2)])
  colnames(serum_fix)[ncol(serum_fix)] <- analytes[k]
}


pdf("hist_test.pdf", 5, 10)
par(mfrow=c(5,3))
for(i in 5:17){
plot(density(serum_fix[,i][!is.na(serum_fix[,i])]), main = colnames(serum_fix)[i])
}
dev.off()

#
#
#
#
#
#
#
#
# check bias

meta <- read.table("../../../../int/data/metadata/mice_metadata.11.19_mouse_id.txt", T, '\t')
cyt_keep <- intersect(serum_fix$sample, meta$mouse_id)

rownames(serum_fix) <- serum_fix$sample
rownames(meta) <- meta$mouse_id

ser_cyt_keep <- serum_fix[cyt_keep,]
ser_cyt_keep[is.na(ser_cyt_keep)] <- 1
ser_cyt_keep[ser_cyt_keep<=1] <- 0
ser_cyt_keep$sample <- rownames(ser_cyt_keep)
meta_keep <- meta[cyt_keep,]

colnames(ser_cyt_keep) <- gsub("Analyte\\.", "", colnames(ser_cyt_keep))
good_cols <- c(5:10,13:16)

write.table(ser_cyt_keep[,c(4,good_cols)],"Serum_names_final.txt", sep='\t', quote = F, row.names=F)

pc <- prcomp(log2(ser_cyt_keep[,good_cols]+1), scale=T)
huh<-pc$rotation

plotter <- data.frame(PC1=pc$x[,1], PC2=pc$x[,2], ser_cyt_keep, meta_keep)

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

table(meta_keep$Gender, ser_cyt_keep$plate)


library(ComplexHeatmap)

pdf("serum_plate_gender_bias.pdf", height = 5, width = 8)
va <- rowAnnotation(df = data.frame(Plate=ser_cyt_keep$plate, Gender=meta_keep$Gender))
Heatmap(log2(ser_cyt_keep[,good_cols]+1), 
        col = colorRamp2(c(1, 7), c("white", "red")),
        cluster_columns = T, cluster_rows = T,
        show_row_names = T, show_column_names = T,
        row_names_gp = gpar(fontsize = 6),
        heatmap_legend_param=list(title = "Log2 Value")) + va
dev.off()


pdf("serum_plate_env.pdf", height = 5, width = 8)
va <- rowAnnotation(df = data.frame(Environment=plotter$Environment, Genotype=plotter$Genotype))
Heatmap(log2(ser_cyt_keep[,good_cols]+1), 
        col = colorRamp2(c(1, 7), c("white", "red")),
        cluster_columns = T, cluster_rows = T,
        show_row_names = T, show_column_names = T,
        row_names_gp = gpar(fontsize = 6),
        heatmap_legend_param=list(title = "Log2 Value")) + va
dev.off()
