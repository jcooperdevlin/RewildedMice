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

## train data X
rnaseq_df <- read.table("normalizedCounts.txt", T, '\t')
rnaseq_var <- rnaseq_df[order(rowVars(data.matrix(rnaseq_df[,-1])),decreasing=T),]#[1:5000,]
trainX <- t(rnaseq_var[,-1])
colnames(trainX) <- rnaseq_var$GeneSymbol
rownames(trainX) <- gsub("X", "", rownames(trainX))
trainXs <- log2(trainX+1)

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

png("gene_cyt_corr_heat.png", height =5, width = 10, units = 'in', res=150)
Heatmap(cor_mat[good_genes,], show_row_names = F, show_column_names = F,
        col = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")))
dev.off()

print(paste0(length(good_genes), " genes selected to train GLM"))
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
  annotate('point', x=best_alpha, y=best_lambda, color="red3", size=3) +
  theme_bw()

pdf("total_loss.pdf", height = 20, width = 10)
gloss
dev.off()

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
  test_x <- trainXs[71:84,]
  test_y <- trainY[71:84,]
  
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

pdf("final_corr.pdf", height = 6, width = 14)
g1
dev.off()

write.table(cor_res, "final_corr_table.txt", sep='\t', quote=F, row.names=F)
save(modelList, file="modelList.RData")


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
###### Now train models with FACS data and maybe 16S?

## train data X
blood_lymph_facs <- read.table("BLOOD_lymph_FACS_metadata_names.txt", T, sep='\t')
remove <- which(duplicated(blood_lymph_facs$mouse_id))
blood_lymph_facs <- blood_lymph_facs[-remove,]
rownames(blood_lymph_facs) <- blood_lymph_facs$mouse_id
colnames(blood_lymph_facs)[13:ncol(blood_lymph_facs)] <- paste0("Blood_", colnames(blood_lymph_facs)[13:ncol(blood_lymph_facs)])

mln_lymph_facs <- read.table("MLN_lymph_FACS_metadata_names.txt", T, sep='\t')
remove <- which(duplicated(mln_lymph_facs$mouse_id))
mln_lymph_facs <- mln_lymph_facs[-remove,]
rownames(mln_lymph_facs) <- mln_lymph_facs$mouse_id
colnames(mln_lymph_facs)[13:ncol(mln_lymph_facs)] <- paste0("MLN_", colnames(mln_lymph_facs)[13:ncol(mln_lymph_facs)])

mln_myeloid_facs <- read.table("MLN_myeloid_FACS_metadata_names.txt", T, sep='\t')
remove <- which(duplicated(mln_myeloid_facs$mouse_id))
mln_myeloid_facs <- mln_myeloid_facs[-remove,]
rownames(mln_myeloid_facs) <- mln_myeloid_facs$mouse_id
colnames(mln_myeloid_facs)[13:ncol(mln_myeloid_facs)] <- paste0("MLN_", colnames(mln_myeloid_facs)[13:ncol(mln_myeloid_facs)])


lymph <- merge(blood_lymph_facs[,c(1,13:ncol(blood_lymph_facs))],
               mln_lymph_facs[,c(1,13:ncol(mln_lymph_facs))],
               by = "mouse_id")
full_facs <- merge(lymph, mln_myeloid_facs[,c(1,13:ncol(mln_myeloid_facs))],
                   by = "mouse_id")
remover <- unique(which(is.na(full_facs), arr.ind = T)[,1])
full_facs <- full_facs[-remover,]
trainXs <- full_facs[,-1]
rownames(trainXs) <- gsub("-", "_", full_facs$mouse_id)
trainXs <- data.matrix(trainXs)


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

png("facs_cyt_corr_heat.png", height =5, width = 10, units = 'in', res=150)
Heatmap(cor_mat[good_genes,], show_row_names = F, show_column_names = F,
        col = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")))
dev.off()

print(paste0(length(good_genes), " FACS channels selected to train GLM"))
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
  annotate('point', x=best_alpha, y=best_lambda, color="red3", size=3) +
  theme_bw()

pdf("total_loss_FACS.pdf", height = 20, width = 10)
gloss
dev.off()

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
  test_x <- trainXs[71:83,]
  test_y <- trainY[71:83,]
  
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

pdf("final_corr_FACS.pdf", height = 6, width = 14)
g1
dev.off()

write.table(cor_res, "final_corr_table_FACS.txt", sep='\t', quote=F, row.names=F)
save(modelList, file="modelList_FACS.RData")

