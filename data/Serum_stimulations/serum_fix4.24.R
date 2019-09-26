#### serum cytokine fix

setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/data/Serum_stimulations")

ser_cyt <- read.table("CBA_Serum.txt", T, '\t')

stim_melt <- melt(ser_cyt, measure.vars = colnames(ser_cyt)[4:16], id.vars = colnames(ser_cyt)[c(1:3,17)])
stim_melt$value <- gsub("<.*", 0, stim_melt$value)
stim_melt$value <- gsub("ND", 0, stim_melt$value)
stim_melt$value <- as.numeric(stim_melt$value)
stim_melt$value[is.na(stim_melt$value)] <- 0
stim_cast <- cast(stim_melt, Name ~ variable, mean)
stim_cast2 <- cast(stim_melt, Name+Genetics+Location+Plate~ variable, mean)

stim_cast$Name <- gsub(" ", "", stim_cast$Name)
stim_cast2$Name <- gsub(" ", "", stim_cast2$Name)

#try to handle bias?

p1 <- subset(stim_cast2, Plate == "P1")
p1_fix <- data.matrix(scale(p1[,5:17]))
p1_fix <- cbind(p1[,1:4], p1_fix)
colnames(p1_fix) <- colnames(p1)
p1_fix[is.na(p1_fix)]<-0

p2 <- subset(stim_cast2, Plate == "P2")
p2_fix <- data.matrix(scale(p2[,5:17]))
p2_fix <- cbind(p2[,1:4], p2_fix)
colnames(p2_fix) <- colnames(p2)
p2_fix[is.na(p2_fix)]<-0

stim_cast_fix <- rbind(p1_fix, p2_fix)
stim_cast <- stim_cast_fix[,c(1,5:17)]

#
#

write.table(stim_cast,"CBA_serum_name.txt", row.names=F, sep='\t', quote=F)



#
#
#
# check bias

meta <- read.table("../../../int/data/metadata/mice_metadata.11.19_mouse_id.txt", T, '\t')
cyt_keep <- intersect(stim_cast_fix$Name, meta$mouse_id)

rownames(stim_cast_fix) <- stim_cast_fix$Name
rownames(meta) <- meta$mouse_id

ser_cyt_keep <- stim_cast_fix[cyt_keep,]
meta_keep <- meta[cyt_keep,]


pc <- prcomp(ser_cyt_keep[,5:17])

plotter <- data.frame(PC1=pc$x[,1], PC2=pc$x[,2], ser_cyt_keep)

ggplot(plotter, aes(PC1, PC2, color = Location)) +
  geom_point()

ggplot(plotter, aes(PC1, PC2, color = Plate)) +
  geom_point()

ggplot(plotter, aes(PC1, PC2, color = Genetics)) +
  geom_point()

ggplot(plotter, aes(PC1, PC2, color = meta_keep$Gender)) +
  geom_point()

table(meta_keep$Gender, ser_cyt_keep$Plate)


library(ComplexHeatmap)

pdf("serum_plate_gender_bias.pdf", height = 5, width = 8)
va <- rowAnnotation(df = data.frame(Plate=ser_cyt_keep$Plate, Gender=meta_keep$Gender))
Heatmap(log2(ser_cyt_keep[,5:17]+1), 
               #col = colorRamp2(c(0, 7), c("white", "red")),
               cluster_columns = T, cluster_rows = T,
               show_row_names = T, show_column_names = T,
               row_names_gp = gpar(fontsize = 6),
               heatmap_legend_param=list(title = "Log2 Value")) + va
dev.off()
