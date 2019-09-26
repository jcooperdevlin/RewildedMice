### Correlation results comparisons

## Best models

setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/IntegrationModel/ModelResults")

library(ggplot2)
library(RColorBrewer)
library(ggsci)

colors_clusters = c(pal_d3("category10")(10), pal_d3("category20b")(20), pal_igv("default")(51))


c_inputs <- list.files(pattern = "final_corr.txt", full.names = T, recursive = T)

c_results <- read.table(c_inputs[1], T, '\t')[1,]
c_results$cytokine=NA
for(i in 1:length(c_inputs)){
  curr <- read.table(c_inputs[i], T, '\t')
  curr$cytokine <- strsplit(c_inputs[i], split = "/")[[1]][2]
  c_results <- rbind(c_results, curr)
  
  
}
c_results <- c_results[-1,]

c_results$condition <- gsub("_Bacillus", ":Bacillus", c_results$condition)
c_results$condition <- gsub("_Bacteroides", ":Bacteroides", c_results$condition)
c_results$condition <- gsub("_Candida", ":Candida", c_results$condition)
c_results$condition <- gsub("_CD3", ":CD3", c_results$condition)
c_results$condition <- gsub("_Clostri", ":Clostri", c_results$condition)
c_results$condition <- gsub("_PBS", ":PBS", c_results$condition)
c_results$condition <- gsub("_Pseudo", ":Pseudo", c_results$condition)
c_results$condition <- gsub("_Staph", ":Staph", c_results$condition)
c_results$Challenge <- gsub(".*:", "", c_results$condition)

g1 <- ggplot(c_results, aes(reorder(Challenge, cor, FUN = median), cor, color=Challenge)) +
  geom_boxplot() + geom_jitter(width=0.1) + ggtitle("Test Correlations") +
  xlab("Condition")+
  scale_color_manual(values=colors_clusters)+
  geom_hline(yintercept = 0.7, color='red') +
  theme_bw() +
  #theme(axis.text.x = element_text(size=10, angle = 90, hjust=1,vjust=1), legend.position = 'right') +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(color='black'))

pdf("Combined_correlation_results_All_models.pdf", height = 10, width = 15)
g1+facet_wrap(~cytokine, ncol=7)
dev.off()


c_good <- aggregate(. ~condition, data=c_results[,1:2], median)
good_good <- c_good$condition[which(c_good$cor >= 0.6)]

c_good_results <- subset(c_results, condition %in% good_good)
g1 <- ggplot(c_good_results, aes(reorder(Challenge, cor, FUN = median), cor, color=Challenge)) +
  geom_boxplot() + geom_jitter(width=0.1) + ggtitle("Test Correlations") +
  xlab("Condition")+
  scale_color_manual(values=colors_clusters)+
  geom_hline(yintercept = 0.7, color='red') +
  theme_bw() +
  #theme(axis.text.x = element_text(size=10, angle = 90, hjust=1,vjust=1), legend.position = 'right') +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(color='black'))

pdf("Combined_correlation_good_results_All_models.pdf", height = 10, width = 15)
g1+facet_wrap(~cytokine, ncol=7)
dev.off()
