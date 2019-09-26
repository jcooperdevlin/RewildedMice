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
  res$color[res$res>0.05] <- "red"
  res$color = factor(res$color, levels = c("black", "red"))
  
  g_res <- ggplot(res, aes(y_act, res, color = color)) + 
    geom_point(show.legend = F) + coord_flip() +
    geom_hline(yintercept = 0.05) +
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
    res$color[res$res>0.05] <- "red"
    res$color = factor(res$color, levels = c("black", "red"))
    
    g_res <- ggplot(res, aes(y_act, res, color = color)) + 
      geom_point(show.legend = F) + coord_flip() +
      geom_hline(yintercept = 0.05) +
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

analytes <- unique(serum_data[[2]]$analyte)
for(k in 1:length(analytes)){
  analyte_name = analytes[k]
  input_mfi <- subset(serum_data[[2]], analyte == analyte_name)
  input_ec <- subset(serum_data[[1]], analyte == analyte_name)
  
  full_plots <- fitter(analyte_name, input_mfi, input_ec)
  namer=paste0("plots/", analyte_name, "_fit-correction.pdf")
  pdf(namer, height = 7, width = 12)
  pp=plot_grid(plotlist = full_plots, ncol=2)
  print(pp)
  dev.off()
}

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
  serum_fix2 <- merge(serum_fix, input_mfi[,c("sample", "est_c")], by="sample", all=F)
  merge(x = DF1, y = DF2[ , c("Client", "LO")], by = "Client", all.x=TRUE)
  
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
# organize data
p1_mfi <- subset(mfi_data, plate=="Plate-1" & analyte == "Analyte.A4.RANTES.")
p1_ec <- subset(ec_data, plate=="Plate-1" & analyte == "Analyte.A4.RANTES.", select=-c(plate))
p1_ec$mfi <- p1_mfi$mfi[1:16]
p1_ec$well <- p1_mfi$well[1:16]

stand_plot <- p1_ec
ggplot(stand_plot, aes(ec+1, mfi+1)) +
  geom_point() + scale_x_log10() + scale_y_log10()

np2 = nplr(log10(stand_plot$ec+1), log10(stand_plot$mfi+1), useLog = F)

x <- getXcurve(np2)
y <-getYcurve(np2)

newdata <- data.frame(x=x,y=y)

ggplot(stand_plot, aes(log10(ec+1), log10(mfi+1))) +
  geom_point() + 
  geom_line(aes(x=x,y=y), newdata)

new_points = getEstimates(np2, log10(stand_plot$mfi+1))

ggplot(stand_plot, aes(log10(ec+1), log10(mfi+1))) +
  geom_point() + 
  geom_line(aes(x=x,y=y), newdata) +
  geom_point(aes(x=x,y=y, col='red'), new_points)

new_points$act_c <- stand_plot$ec
new_points$est_c <- 10^new_points$x
new_points$recovery <- (new_points$act_c/new_points$est_c)*100

myPalette <- colorRampPalette("dodgerblue2")
ggplot(stand_plot, aes(log10(ec+1), log10(mfi+1))) +
  geom_point() + 
  geom_line(aes(x=x,y=y), newdata) +
  geom_point(aes(x=x,y=y, col=recovery), new_points) +
  scale_colour_gradientn(colours = myPalette(n=10), limits=c(80, 120))

##
#
#
#
### remove some points
p1_mfi <- subset(mfi_data, plate=="Plate-1" & analyte == "Analyte.A4.RANTES.")
p1_ec <- subset(ec_data, plate=="Plate-1" & analyte == "Analyte.A4.RANTES.", select=-c(plate))
p1_ec$mfi <- p1_mfi$mfi[1:16]
p1_ec$well <- p1_mfi$well[1:16]

p1_mfi <- subset(p1_mfi, sample != "Standard3")
p1_ec <- subset(p1_ec, sample != "Standard3")

stand_plot <- p1_ec
ggplot(stand_plot, aes(ec+1, mfi+1)) +
  geom_point() + scale_x_log10() + scale_y_log10()

np2 = nplr(log10(stand_plot$ec+1), log10(stand_plot$mfi+1), useLog = F)

x <- getXcurve(np2)
y <-getYcurve(np2)

newdata <- data.frame(x=x,y=y)

ggplot(stand_plot, aes(log10(ec+1), log10(mfi+1))) +
  geom_point() + 
  geom_line(aes(x=x,y=y), newdata)

new_points = getEstimates(np2, log10(stand_plot$mfi+1))

ggplot(stand_plot, aes(log10(ec+1), log10(mfi+1))) +
  geom_point() + 
  geom_line(aes(x=x,y=y), newdata) +
  geom_point(aes(x=x,y=y, col='red'), new_points)

new_points$act_c <- stand_plot$ec
new_points$est_c <- 10^new_points$x
new_points$recovery <- (new_points$act_c/new_points$est_c)*100

myPalette <- colorRampPalette("dodgerblue2")
ggplot(stand_plot, aes(log10(ec+1), log10(mfi+1))) +
  geom_point() + 
  geom_line(aes(x=x,y=y), newdata) +
  geom_point(aes(x=x,y=y, col=recovery), new_points) +
  scale_colour_gradientn(colours = myPalette(n=10), limits=c(80, 120))


#
#
### remove some points and avg replicates?

### remove some points
p1_mfi <- subset(mfi_data, plate=="Plate-1" & analyte == "Analyte.A4.RANTES.")
p1_mfi <- subset(p1_mfi, grepl("Background", p1_mfi$sample)==F)
p1_ec <- subset(ec_data, plate=="Plate-1" & analyte == "Analyte.A4.RANTES.", select=-c(plate))
p1_ec$mfi <- p1_mfi$mfi[1:16]

p1_mfi <- subset(p1_mfi, sample != "Standard3")
p1_ec <- subset(p1_ec, sample != "Standard3")

p1_ec_agg <- aggregate(mfi ~ sample + analyte + ec, mean, data=p1_ec)

stand_plot <- p1_ec_agg
ggplot(stand_plot, aes(ec+1, mfi+1)) +
  geom_point() + scale_x_log10() + scale_y_log10()

np2 = nplr(log10(stand_plot$ec+1), log10(stand_plot$mfi+1), useLog = F)
plot(np2)

x <- getXcurve(np2)
y <-getYcurve(np2)

newdata <- data.frame(x=x,y=y)

ggplot(stand_plot, aes(log10(ec+1), log10(mfi+1))) +
  geom_point() + 
  geom_line(aes(x=x,y=y), newdata)

new_points = getEstimates(np2, log10(stand_plot$mfi+1))

ggplot(stand_plot, aes(log10(ec+1), log10(mfi+1))) +
  geom_point() + 
  geom_line(aes(x=x,y=y), newdata) +
  geom_point(aes(x=x,y=y, col='red'), new_points)

new_points$act_c <- stand_plot$ec
new_points$est_c <- 10^new_points$x
new_points$recovery <- (new_points$act_c/new_points$est_c)*100

myPalette <- colorRampPalette("dodgerblue2")
ggplot(stand_plot, aes(log10(ec+1), log10(mfi+1))) +
  geom_point() + 
  geom_line(aes(x=x,y=y), newdata) +
  geom_point(aes(x=x,y=y, col=recovery), new_points) +
  scale_colour_gradientn(colours = myPalette(n=10), limits=c(80, 120))

#
#
#

#
### estimate the unknowns based on chosen fit
p1_mfi$mfi_log <- log10(p1_mfi$mfi+1)
est_conc <- getEstimates(np2, p1_mfi$mfi_log)
p1_mfi$est_ec_log <-est_conc$x
p1_mfi$est_ec <-10^est_conc$x

ggplot(stand_plot, aes(log10(ec+1), log10(mfi+1))) +
  geom_point() + 
  geom_line(aes(x=x,y=y), newdata) +
  geom_point(aes(x=x,y=y, col=recovery), new_points) +
  scale_colour_gradientn(colours = myPalette(n=10), limits=c(80, 120)) +
  geom_point(aes(x=est_ec_log, y = mfi_log), p1_mfi)
##


p1_mfi <- subset(p1_mfi, sample != "Standard3")
p1_ec <- subset(p1_ec, sample != "Standard3")

p1_ec2 <- p1_ec[c(1,3,5,7,9,11,13,15),]


stand_plot <- p1_ec2
ggplot(stand_plot, aes(ec+1, mfi+1)) +
  geom_point() + scale_x_log10() + scale_y_log10()




fit <- scluminex(plateid="Plate-1", standard = p1_ec2, 
                 background = data.frame(),
                 bkg = "ignore", lfct = c("SSl5","SSl4"),
                 fmfi = "mfi", verbose = TRUE)


datasets <- data_selection(x = p1_mfi, ecfile = p1_ec,
                           byvar.ecfile = c("sample","analyte"),
                           backname = "Background0",
                           stanname="Standard")

stand_plot <- datasets$`Plate-1`$standard
ggplot(stand_plot, aes(ec+1, mfi+1)) +
  geom_point() + scale_x_log10() + scale_y_log10()

p1_mfi <- subset(p1_mfi, sample != "Standard5" & sample != "Standard6")
p1_ec <- subset(p1_ec, sample != "Standard5" & sample != "Standard6")

datasets <- data_selection(x = p1_mfi, ecfile = p1_ec,
                           byvar.ecfile = c("sample","analyte"),
                           backname = "Background0",
                           stanname="Standard")

stand_plot <- datasets$`Plate-1`$standard
ggplot(stand_plot, aes(ec+1, mfi+1)) +
  geom_point() + scale_x_log10() + scale_y_log10()

allanalytes <- scluminex(plateid = "newplate",
                         standard = datasets$`Plate-1`$standard,
                         background = datasets$`Plate-1`$background,
                         bkg = "ignore", lfct = c("SSl5","SSl4"),
                         fmfi = "mfi", verbose = TRUE)





summary(allanalytes)



#
#
#
#tests
data(mfidata)
data(ecdata)

mfidata$well = NA
mfidata <- subset(mfidata, grepl("Control", mfidata$sample)==F)
ecdata <- subset(ecdata, grepl("Control", ecdata$sample)==F)

ec_simp <- subset(ecdata, analyte == "FGF")
mfi_simp <- subset(mfidata, analyte == "FGF" & grepl("Standard", mfidata$sample))

ec_simp=ec_simp[order(ec_simp$sample, decreasing=T),]
mfi_simp=mfi_simp[order(mfi_simp$sample, decreasing=T),]

ec_simp$mfi <- c(mfi_simp$mfi, 0,0,0)
ec_simp <- ec_simp[1:16,]
ec_simp$well = NA


fit <- scluminex(plateid="newplate", standard = ec_simp, 
                 background = data.frame(),
                 bkg = "ignore", lfct = c("SSl5","SSl4"),
                 fmfi = "mfi", verbose = TRUE)



datasets <- data_selection(x = mfidata, ecfile = ecdata,
                           byvar.ecfile = c("sample","analyte"),
                           #backname = "Background0",
                           stanname="Standard")



allanalytes <- scluminex(plateid = "newplate",
                         standard = datasets$plate_1$standard,
                         background = datasets$plate_1$background,
                         bkg = "ignore", lfct = c("SSl5","SSl4"),
                         fmfi = "mfi", verbose = FALSE)




#
#
#

#
#
#
data(ecdata)
data(mfidata)

dat <- mfidata[mfidata$plate=="plate_1" & mfidata$analyte=="FGF",]

sdf <- data_selection(dat, ecdata)$plate_1

# Fit model and summary object
igmodels <- scluminex("plate_1",sdf$standard, sdf$background,
                      lfct=c("SSl4", "SSl5"),
                      bkg="ignore",
                      fmfi="mfi",
                      verbose=FALSE)
ss <- summary(igmodels)

# Information
names(igmodels)
names(igmodels$FGF)

# Summary data
ss
as.data.frame(ss)
as.data.frame(igmodels)

# Plot the standard curve
plot(igmodels,"sc")


#
#
#
#
#
#
#extra stuff
sc <- data.frame(Name=NA, Plate=NA, Concentration=NA, MFI=NA, c_log=NA, MFI_log=NA)
for (i in 4:16){
  name <- colnames(read.xlsx("SerumCBA-Plate1.detail.xlsx", sheetIndex = i, startRow = 0, endRow = 22))[1]
  x <- read.xlsx("SerumCBA-Plate1.detail.xlsx", sheetIndex = i, startRow = 5, endRow = 22)
  x$Plate <- rep("Plate-1", nrow(x))
  x$Name <- rep(name, nrow(x))
  exp <- unique(x$Expected.pg.ml.)
  x$Concentration <- rep(exp[c(1,3:9)], each = 2)
  
  x$c_log <- log10(x$Concentration+1)
  x$MFI_log <- log10(x$MFI+1)
  
  x_add <- x[,c("Name", "Plate", "Concentration", "MFI", "c_log", "MFI_log")]
  
  sc <- rbind(sc, x_add)
  
  name <- colnames(read.xlsx("SerumCBA-Plate1.detail.xlsx", sheetIndex = i, startRow = 0, endRow = 22))[1]
  x <- read.xlsx("SerumCBA-Plate2.detail.xlsx", sheetIndex = i, startRow = 5, endRow = 22)
  x$Plate <- rep("Plate-2", nrow(x))
  x$Name <- rep(name, nrow(x))
  exp <- unique(x$Expected.pg.ml.)
  x$Concentration <- rep(exp[c(1,3:9)], each = 2)
  
  x$c_log <- log10(x$Concentration+1)
  x$MFI_log <- log10(x$MFI+1)
  
  x_add <- x[,c("Name", "Plate", "Concentration", "MFI", "c_log", "MFI_log")]
  
  sc <- rbind(sc, x_add)
}
sc <- sc[-1,]

cyts <- unique(sc$Name)

for(i in 1:length(cyts)){
  
  
  p1 <- subset(sc, Plate == "Plate-1" & Name == cyts[i])
  p2 <- subset(sc, Plate == "Plate-2" & Name == cyts[i])
  
  
  yp <- convertToProp(p1$Concentration)
  np1 <- nplr(x=p1$MFI, y=yp, useLog = T)
  
  plot(np1)
  
  
  
  fit.5pl.p1 <- drm(c_log ~ MFI_log, data = p1, fct = LL.5())
  summary(fit.5pl.p1)
  plot(fit.5pl.p1, col=rep(c("red"), nrow(p1)), add=F)
  
  fit.5pl.p2 <- drm(c_log ~ MFI_log, data = p2, fct = LL.5())
  summary(fit.5pl.p2)
  plot(fit.5pl.p2, col=rep(c("blue"), nrow(p2)), add=T)
  
  p12 <- rbind(p1,p2)
  fit.5pl.p12 <- drm(c_log ~ MFI_log, data = p12, fct = LL.5())
  summary(fit.5pl.p12)
  plot(fit.5pl.p12, col=rep(c("purple"), nrow(p12)), add=F)
  
  p12$model_guess <- predict(fit.5pl.p12, data.frame(p12$MFI_log))
  p12$model_conc <- 10^p12$model_guess
  
  
  p1 <- subset(sc, Plate == "Plate-1" & Name == cyts[i] & Concentration != "9.77")
  p2 <- subset(sc, Plate == "Plate-2" & Name == cyts[i] & Concentration != "9.77")
  
  p12 <- rbind(p1,p2)
  fit.5pl.p12 <- drm(c_log ~ MFI_log, data = p12, fct = LL.5())
  summary(fit.5pl.p12)
  plot(fit.5pl.p12, col=rep(c("purple"), nrow(p12)), add=F)
  
  p12$model_guess <- predict(fit.5pl.p12, data.frame(p12$MFI_log))
  p12$model_conc <- 10^p12$model_guess
  
}
