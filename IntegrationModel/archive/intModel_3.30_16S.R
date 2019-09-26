### Integration Model Test1
### 3.27

#setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/IntegrationModel")

library(ggplot2)
#library(MultivariateRandomForest)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(psych)
library(gridExtra)

source("glmnet_full_function.R")
###### Now train models with 16S

## train data X
otu_tab <- read.table("otu_table_counts.txt", T, '\t')
#otu_tab <- aggregate(otu_tab, by=list("taxonomy"), sum)
otu_tab <- aggregate(.~taxonomy, data=otu_tab, sum)
rownames(otu_tab) <- otu_tab$taxonomy
otu_tab <- otu_tab[,-1]
otu_tab <- t(otu_tab)
otu_tab <- otu_tab[,-1]

rownames(otu_tab) <- gsub("X", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("_A.*", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("_N.*", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("_B.*", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("_", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("\\.", "_", rownames(otu_tab))

trainXs <- log2(otu_tab+1)


### train data y
cyt_df <- read.table("MLN_stimulation_flat.txt", T, '\t')
rownames(cyt_df) <- cyt_df$mouse_id
cyt_df <- log2(cyt_df[,-1]+1)
#cyt_df <- scale(cyt_df[,-1])
trainY <- data.matrix(cyt_df)
rownames(trainY) <- gsub("-", "_", rownames(trainY))

overlaps <- intersect(rownames(trainY), rownames(trainXs))
overlaps <- sample(overlaps, length(overlaps))

trainXs <- trainXs[overlaps,]
trainY <- trainY[overlaps,]


run_glm(trainXs, trainY, "test", type="OTUs", p_threshold = 0.05)


### are predictors correlated?

cor_adder <- c()
p_adder <- c()
for(i in 1:ncol(trainY)){
  cc <- apply(trainXs[overlaps,], 2, function(x) {cor.test(x, trainY[overlaps,i], ci = F)})
  rr <- round(unlist(sapply(cc, `[`, 4)),3)
  pp <- round(unlist(sapply(cc, `[`, 3)),3)
  cor_adder <- c(cor_adder, rr)
  p_adder <- c(p_adder, pp)
  if(i %% 10==0){print(paste0(i, " stimulations tested out of ", ncol(trainY)))}
}

cor_mat <- matrix(cor_adder, ncol(trainXs), 104)
rownames(cor_mat) <- colnames(trainXs)
colnames(cor_mat) <- colnames(trainY)

p_adder <- p.adjust(p_adder, method="BH")
p_mat <- matrix(p_adder, ncol(trainXs), 104)
rownames(p_mat) <- colnames(trainXs)
colnames(p_mat) <- colnames(trainY)

good_genes <- rownames(cor_mat)[unique(which(p_mat <= 0.05, arr.ind = T)[,1])]

png("otu_cyt_corr_heat.png", height =5, width = 10, units = 'in', res=150)
Heatmap(cor_mat[good_genes,], show_row_names = F, show_column_names = F,
        col = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")))
dev.off()

print(paste0(length(good_genes), " OTUs selected to train GLM"))
#
#

### RANDOMLY TEST 10 Times different test and train sets for full dataset with all genes

library(caret)
library(glmnet)

results_df <- data.frame(alpha=NA, lambda=NA,MSE=NA, num=NA)
alpha_vals <- seq(0,1, length=50)
#lambda_vals <- seq(0.001,0.1,by = 0.001)
lambda_vals <- seq(0.0001, 1, length = 100)
trainModelList=list()
for (k in 1:50) {
  if(k %% 5==0){print(paste0("Starting iteration...", k))}
  trainModelList[[k]]<- cv.glmnet(trainXs, trainY, type.measure="mse",nfolds = 7,
                                  alpha=alpha_vals[k],family="mgaussian", #lambda = lambda_vals,
                                  standardize.response=T)
  
  adder <- data.frame(alpha = rep(round(alpha_vals[k],3),length(trainModelList[[k]]$lambda)), lambda = trainModelList[[k]]$lambda, 
                      MSE = trainModelList[[k]]$cvm, num=k)
  results_df <- rbind(results_df, adder)
}
results_df<-results_df[-1,]
results_df <- results_df[order(results_df$MSE, decreasing=F),]
bestMod <- results_df$num[1]

bestModel <- trainModelList[[bestMod]]
best_alpha <- results_df$alpha[1]
best_lambda <- results_df$lambda[1]

results_df$Log2_MSE <- log2(results_df$MSE+1)
gloss <- ggplot(results_df, aes(x=alpha, y=log2(lambda+1), color = Log2_MSE))+
  geom_jitter(width=0.01,size=1) + ylab("Log2 lambda") +
  scale_color_gradient2(low="purple4", mid="mediumorchid3", high = "grey80", 
                        midpoint=median(results_df$Log2_MSE)) +
  annotate('point', x=best_alpha, y=log2(best_lambda+1), color="red3", size=3) +
  theme_bw()

pdf("total_loss_16S.pdf", height = 10, width = 10)
gloss
dev.off()

write.table(results_df, "loss_results_16S.txt", row.names=F, sep='\t', quote=F)

print(paste0("Alpha selected at...", best_alpha))
print(paste0("Lambda selected at...", round(best_lambda,3)))

cor_res <- data.frame(condition=NA, cor=NA, iteration=NA)
modelList <- list()
for(i in 1:10){
  ## check intersect
  print(paste0("Starting CV iteration...", i))
  
  overlaps <- intersect(rownames(trainY), rownames(trainXs))
  overlaps <- sample(overlaps, length(overlaps))
  
  trainXs <- trainXs[overlaps,good_genes]
  trainY <- trainY[overlaps,]
  
  train_x <- trainXs[1:70,]
  train_y <- trainY[1:70,]
  test_x <- trainXs[71:88,]
  test_y <- trainY[71:88,]
  
  modelList[[i]]<- glmnet(train_x, train_y, #type.measure="mse",nfolds = 7,
                          alpha=best_alpha,family="mgaussian", lambda = best_lambda,
                          standardize.response=T)
  
  yhat <- matrix(predict(modelList[[i]], s=modelList[[i]]$lambda, newx=test_x), 
                 nrow(test_x),ncol(test_y))
  mse <- mean((test_y - yhat)^2)
  
  for (j in 1:ncol(test_y)){
    condition<-colnames(test_y)[j]
    cc=cor(yhat[,j], test_y[,j])
    iteration=i
    adder <- data.frame(condition = condition, cor = cc, iteration = iteration)
    cor_res <- rbind(cor_res, adder)
  }
}
cor_res <- cor_res[-1,]

g1 <- ggplot(cor_res, aes(reorder(condition, cor, FUN = median), cor, color=condition)) +
  geom_boxplot() + geom_jitter(width=0.1) + ggtitle("Test Correlations") +
  xlab("Condition")+
  theme(axis.text.x = element_text(size=10, angle = 90, hjust=1,vjust=1), legend.position = 'none')

pdf("final_corr_OTUs.pdf", height = 6, width = 14)
g1
dev.off()

write.table(cor_res, "final_corr_table_OTUs.txt", sep='\t', quote=F, row.names=F)
save(modelList, file="modelList_OTUs.RData")


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
###### Now train models with 16S

## train data X
otu_tab <- read.table("otu_table.txt", T, '\t')
otu_tab <- aggregate(otu_tab, by=list("taxonomy"), sum)
otu_tab <- aggregate(.~taxonomy, data=otu_tab, sum)
rownames(otu_tab) <- otu_tab$taxonomy
otu_tab <- otu_tab[,-1]
otu_tab <- t(otu_tab)
otu_tab <- otu_tab[,-1]

rownames(otu_tab) <- gsub("X", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("_A.*", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("_N.*", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("_B.*", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("_", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("\\.", "_", rownames(otu_tab))

trainXs <- otu_tab


### train data y
cyt_df <- read.table("MLN_stimulation_flat.txt", T, '\t')
rownames(cyt_df) <- cyt_df$mouse_id
cyt_df <- log2(cyt_df[,-1]+1)
#cyt_df <- scale(cyt_df[,-1])
trainY <- data.matrix(cyt_df)
rownames(trainY) <- gsub("-", "_", rownames(trainY))

overlaps <- intersect(rownames(trainY), rownames(trainXs))
overlaps <- sample(overlaps, length(overlaps))

trainXs <- trainXs[overlaps,]
trainY <- trainY[overlaps,]

### are predictors correlated?

cor_adder <- c()
p_adder <- c()
for(i in 1:ncol(trainY)){
  cc <- apply(trainXs[overlaps,], 2, function(x) {cor.test(x, trainY[overlaps,i], ci = F)})
  rr <- round(unlist(sapply(cc, `[`, 4)),3)
  pp <- round(unlist(sapply(cc, `[`, 3)),3)
  cor_adder <- c(cor_adder, rr)
  p_adder <- c(p_adder, pp)
  if(i %% 10==0){print(paste0(i, " stimulations tested out of ", ncol(trainY)))}
}

cor_mat <- matrix(cor_adder, ncol(trainXs), 104)
rownames(cor_mat) <- colnames(trainXs)
colnames(cor_mat) <- colnames(trainY)

p_adder <- p.adjust(p_adder, method="BH")
p_mat <- matrix(p_adder, ncol(trainXs), 104)
rownames(p_mat) <- colnames(trainXs)
colnames(p_mat) <- colnames(trainY)

good_genes <- rownames(cor_mat)[unique(which(p_mat <= 0.05, arr.ind = T)[,1])]

png("otu_cyt_corr_heat.png", height =5, width = 10, units = 'in', res=150)
Heatmap(cor_mat[good_genes,], show_row_names = F, show_column_names = F,
        col = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")))
dev.off()

print(paste0(length(good_genes), " OTUs selected to train GLM"))
#
#

### RANDOMLY TEST 10 Times different test and train sets for full dataset with all genes

library(caret)
library(glmnet)

results_df <- data.frame(alpha=NA, lambda=NA,MSE=NA, num=NA)
alpha_vals <- seq(0,1, length=50)
#lambda_vals <- seq(0.001,0.1,by = 0.001)
lambda_vals <- seq(0.0001, 1, length = 100)
trainModelList=list()
for (k in 1:50) {
  if(k %% 5==0){print(paste0("Starting iteration...", k))}
  trainModelList[[k]]<- cv.glmnet(trainXs, trainY, type.measure="mse",nfolds = 7,
                                  alpha=alpha_vals[k],family="mgaussian", #lambda = lambda_vals,
                                  standardize.response=T)
  
  adder <- data.frame(alpha = rep(round(alpha_vals[k],3),length(trainModelList[[k]]$lambda)), lambda = trainModelList[[k]]$lambda, 
                      MSE = trainModelList[[k]]$cvm, num=k)
  results_df <- rbind(results_df, adder)
}
results_df<-results_df[-1,]
results_df <- results_df[order(results_df$MSE, decreasing=F),]
bestMod <- results_df$num[1]

bestModel <- trainModelList[[bestMod]]
best_alpha <- results_df$alpha[1]
best_lambda <- results_df$lambda[1]

results_df$Log2_MSE <- log2(results_df$MSE+1)
gloss <- ggplot(results_df, aes(x=alpha, y=log2(lambda+1), color = Log2_MSE))+
  geom_jitter(width=0.01,size=1) + ylab("Log2 lambda") +
  scale_color_gradient2(low="purple4", mid="mediumorchid3", high = "grey80", 
                        midpoint=median(results_df$Log2_MSE)) +
  annotate('point', x=best_alpha, y=log2(best_lambda+1), color="red3", size=3) +
  theme_bw()

pdf("total_loss_16S.pdf", height = 10, width = 10)
gloss
dev.off()

write.table(results_df, "loss_results_16S.txt", row.names=F, sep='\t', quote=F)

print(paste0("Alpha selected at...", best_alpha))
print(paste0("Lambda selected at...", round(best_lambda,3)))

cor_res <- data.frame(condition=NA, cor=NA, iteration=NA)
modelList <- list()
for(i in 1:10){
  ## check intersect
  print(paste0("Starting CV iteration...", i))
  
  overlaps <- intersect(rownames(trainY), rownames(trainXs))
  overlaps <- sample(overlaps, length(overlaps))
  
  trainXs <- trainXs[overlaps,good_genes]
  trainY <- trainY[overlaps,]
  
  train_x <- trainXs[1:70,]
  train_y <- trainY[1:70,]
  test_x <- trainXs[71:88,]
  test_y <- trainY[71:88,]
  
  modelList[[i]]<- glmnet(train_x, train_y, #type.measure="mse",nfolds = 7,
                          alpha=best_alpha,family="mgaussian", lambda = best_lambda,
                          standardize.response=T)
  
  yhat <- matrix(predict(modelList[[i]], s=modelList[[i]]$lambda, newx=test_x), 
                 nrow(test_x),ncol(test_y))
  mse <- mean((test_y - yhat)^2)
  
  for (j in 1:ncol(test_y)){
    condition<-colnames(test_y)[j]
    cc=cor(yhat[,j], test_y[,j])
    iteration=i
    adder <- data.frame(condition = condition, cor = cc, iteration = iteration)
    cor_res <- rbind(cor_res, adder)
  }
}
cor_res <- cor_res[-1,]

g1 <- ggplot(cor_res, aes(reorder(condition, cor, FUN = median), cor, color=condition)) +
  geom_boxplot() + geom_jitter(width=0.1) + ggtitle("Test Correlations") +
  xlab("Condition")+
  theme(axis.text.x = element_text(size=10, angle = 90, hjust=1,vjust=1), legend.position = 'none')

pdf("final_corr_OTUs.pdf", height = 6, width = 14)
g1
dev.off()

write.table(cor_res, "final_corr_table_OTUs.txt", sep='\t', quote=F, row.names=F)
save(modelList, file="modelList_OTUs.RData")
#

##
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
#
#

### relab
###### Now train models with 16S

## train data X
otu_tab <- read.table("otu_table.txt", T, '\t', skip=1)
rownames(otu_tab) <- otu_tab$ID
otu_tab <- otu_tab[,-1]
otu_tab <- t(otu_tab)
otu_tab <- otu_tab[,-1]

rownames(otu_tab) <- gsub("X", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("_A.*", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("_N.*", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("_B.*", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("_", "", rownames(otu_tab))
rownames(otu_tab) <- gsub("\\.", "_", rownames(otu_tab))

trainXs <- otu_tab


### train data y
cyt_df <- read.table("MLN_stimulation_flat.txt", T, '\t')
rownames(cyt_df) <- cyt_df$mouse_id
cyt_df <- log2(cyt_df[,-1]+1)
#cyt_df <- scale(cyt_df[,-1])
trainY <- data.matrix(cyt_df)
rownames(trainY) <- gsub("-", "_", rownames(trainY))

overlaps <- intersect(rownames(trainY), rownames(trainXs))
overlaps <- sample(overlaps, length(overlaps))

trainXs <- trainXs[overlaps,]
trainY <- trainY[overlaps,]

### are predictors correlated?

cor_adder <- c()
p_adder <- c()
for(i in 1:ncol(trainY)){
  cc <- apply(trainXs[overlaps,], 2, function(x) {cor.test(x, trainY[overlaps,i], ci = F)})
  rr <- round(unlist(sapply(cc, `[`, 4)),3)
  pp <- round(unlist(sapply(cc, `[`, 3)),3)
  cor_adder <- c(cor_adder, rr)
  p_adder <- c(p_adder, pp)
  if(i %% 10==0){print(paste0(i, " stimulations tested out of ", ncol(trainY)))}
}

cor_mat <- matrix(cor_adder, ncol(trainXs), 104)
rownames(cor_mat) <- colnames(trainXs)
colnames(cor_mat) <- colnames(trainY)

p_adder <- p.adjust(p_adder, method="BH")
p_mat <- matrix(p_adder, ncol(trainXs), 104)
rownames(p_mat) <- colnames(trainXs)
colnames(p_mat) <- colnames(trainY)

good_genes <- rownames(cor_mat)[unique(which(p_mat <= 0.05, arr.ind = T)[,1])]

png("otu_cyt_corr_heat.png", height =5, width = 10, units = 'in', res=150)
Heatmap(cor_mat[good_genes,], show_row_names = F, show_column_names = F,
        col = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")))
dev.off()

print(paste0(length(good_genes), " OTUs selected to train GLM"))
#
#

### RANDOMLY TEST 10 Times different test and train sets for full dataset with all genes

library(caret)
library(glmnet)

results_df <- data.frame(alpha=NA, lambda=NA,MSE=NA, num=NA)
alpha_vals <- seq(0,1, length=50)
#lambda_vals <- seq(0.001,0.1,by = 0.001)
lambda_vals <- seq(0.0001, 1, length = 100)
trainModelList=list()
for (k in 1:50) {
  if(k %% 5==0){print(paste0("Starting iteration...", k))}
  trainModelList[[k]]<- cv.glmnet(trainXs, trainY, type.measure="mse",nfolds = 7,
                                  alpha=alpha_vals[k],family="mgaussian", #lambda = lambda_vals,
                                  standardize.response=T)
  
  adder <- data.frame(alpha = rep(round(alpha_vals[k],3),length(trainModelList[[k]]$lambda)), lambda = trainModelList[[k]]$lambda, 
                      MSE = trainModelList[[k]]$cvm, num=k)
  results_df <- rbind(results_df, adder)
}
results_df<-results_df[-1,]
results_df <- results_df[order(results_df$MSE, decreasing=F),]
bestMod <- results_df$num[1]

bestModel <- trainModelList[[bestMod]]
best_alpha <- results_df$alpha[1]
best_lambda <- results_df$lambda[1]

results_df$Log2_MSE <- log2(results_df$MSE+1)
gloss <- ggplot(results_df, aes(x=alpha, y=log2(lambda+1), color = Log2_MSE))+
  geom_jitter(width=0.01,size=1) + ylab("Log2 lambda") +
  scale_color_gradient2(low="purple4", mid="mediumorchid3", high = "grey80", 
                        midpoint=median(results_df$Log2_MSE)) +
  annotate('point', x=best_alpha, y=log2(best_lambda+1), color="red3", size=3) +
  theme_bw()

pdf("total_loss_FACS.pdf", height = 20, width = 10)
gloss
dev.off()

write.table(results_df, "loss_results_16S.txt", row.names=F, sep='\t', quote=F)

print(paste0("Alpha selected at...", best_alpha))
print(paste0("Lambda selected at...", round(best_lambda,3)))

cor_res <- data.frame(condition=NA, cor=NA, iteration=NA)
modelList <- list()
for(i in 1:10){
  ## check intersect
  print(paste0("Starting CV iteration...", i))
  
  overlaps <- intersect(rownames(trainY), rownames(trainXs))
  overlaps <- sample(overlaps, length(overlaps))
  
  trainXs <- trainXs[overlaps,good_genes]
  trainY <- trainY[overlaps,]
  
  train_x <- trainXs[1:70,]
  train_y <- trainY[1:70,]
  test_x <- trainXs[71:88,]
  test_y <- trainY[71:88,]
  
  modelList[[i]]<- glmnet(train_x, train_y, #type.measure="mse",nfolds = 7,
                          alpha=best_alpha,family="mgaussian", lambda = best_lambda,
                          standardize.response=T)
  
  yhat <- matrix(predict(modelList[[i]], s=modelList[[i]]$lambda, newx=test_x), 
                 nrow(test_x),ncol(test_y))
  mse <- mean((test_y - yhat)^2)
  
  for (j in 1:ncol(test_y)){
    condition<-colnames(test_y)[j]
    cc=cor(yhat[,j], test_y[,j])
    iteration=i
    adder <- data.frame(condition = condition, cor = cc, iteration = iteration)
    cor_res <- rbind(cor_res, adder)
  }
}
cor_res <- cor_res[-1,]

g1 <- ggplot(cor_res, aes(reorder(condition, cor, FUN = median), cor, color=condition)) +
  geom_boxplot() + geom_jitter(width=0.1) + ggtitle("Test Correlations") +
  xlab("Condition")+
  theme(axis.text.x = element_text(size=10, angle = 90, hjust=1,vjust=1), legend.position = 'none')

pdf("final_corr_OTUs.pdf", height = 6, width = 14)
g1
dev.off()

write.table(cor_res, "final_corr_table_OTUs.txt", sep='\t', quote=F, row.names=F)
save(modelList, file="modelList_OTUs.RData")

