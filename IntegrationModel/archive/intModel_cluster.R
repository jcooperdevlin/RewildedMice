### Integration Model Test1
### 3.27

setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/IntegrationModel")

library(ggplot2)
#library(MultivariateRandomForest)
library(matrixStats)
library(ComplexHeatmap)
library(gridExtra)
library(circlize)

## train data X
rnaseq_df <- read.table("normalizedCounts.txt", T, '\t')
rnaseq_var <- rnaseq_df[order(rowVars(data.matrix(rnaseq_df[,-1])),decreasing=T),]#[1:2000,]
trainX <- t(rnaseq_var[,-1])
colnames(trainX) <- rnaseq_var$GeneSymbol
rownames(trainX) <- gsub("X", "", rownames(trainX))
trainXs <- log2(trainX+1)

### train data y
cyt_df <- read.table("MLN_stimulation_flat.txt", T, '\t')
rownames(cyt_df) <- cyt_df$mouse_id
cyt_df <- log2(cyt_df[,-1]+1)
trainY <- data.matrix(cyt_df)
rownames(trainY) <- gsub("-", "_", rownames(trainY))


##
overlaps <- intersect(rownames(trainY), rownames(trainXs))
overlaps <- sample(overlaps, length(overlaps))


cor_adder <- c()
for(i in 1:ncol(trainY)){
  cc <- apply(trainXs[overlaps,], 2, function(x) {cor(x, trainY[overlaps,i])})
  cor_adder <- c(cor_adder, cc)
  if(i %% 10==0){print(paste0(i, " stimulations tested out of ", ncol(trainY)))}
}

cor_mat <- matrix(cor_adder, 104, ncol(trainXs))
colnames(cor_mat) <- colnames(trainXs)
rownames(cor_mat) <- colnames(trainY)

pdf("gene_cyt_corr_heat.pdf", height =5, width = 8)
Heatmap(cor_mat, show_row_names = F, show_column_names = F,
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))
dev.off()

good_genes <- colnames(cor_mat)[unique(which(cor_mat >= 0.5, arr.ind = T)[,2])]

print(paste0(length(good_genes), " genes selected to train GLM"))
#
#

### RANDOMLY TEST 10 Times different test and train sets for full dataset with all genes
cor_res_train <- data.frame(condition=NA, cor=NA, iteration=NA)
cor_res_test <- data.frame(condition=NA, cor=NA, iteration=NA)
modelList <- list()
results_df <- data.frame(alpha=NA, lambda=NA,MSE=NA, num=NA, iteration=NA)
for(i in 1:10){
  ## check intersect
  print(paste0("Starting iteration...", i))
  
  overlaps <- intersect(rownames(trainY), rownames(trainXs))
  overlaps <- sample(overlaps, length(overlaps))
  
  trainXs <- trainXs[overlaps,good_genes]
  trainY <- trainY[overlaps,]
  
  train_x <- trainXs[1:70,]
  train_y <- trainY[1:70,]
  test_x <- trainXs[71:84,]
  test_y <- trainY[71:84,]
  
  library(caret)
  library(glmnet)
  
  # 10-fold Cross validation for each alpha = 0, 0.1, ... , 0.9, 1.0
  # (For plots on Right)
  temp_results_df <- data.frame(alpha=NA, lambda=NA,MSE=NA, num=NA, iteration=NA)
  alpha_vals <- seq(0,1, length=50)
  trainModelList=list()
  for (k in 1:50) {
    if (k %% 10==0){print(paste0("Working on alpha value...", k, " of ", length(alpha_vals)))}
    trainModelList[[k]]<- cv.glmnet(train_x, train_y, type.measure="mse",nfolds = 7,
                                    alpha=alpha_vals[k],family="mgaussian", 
                                    standardize.response=T)
    yhat <- matrix(predict(trainModelList[[k]], s=trainModelList[[k]]$lambda.1se, newx=test_x), 
                   nrow(test_x),ncol(test_y))
    mse <- mean((test_y - yhat)^2)
    adder <- data.frame(alpha = round(alpha_vals[k],3), lambda = trainModelList[[k]]$lambda.1se, 
                        MSE = mse, num=k, iteration=i)
    temp_results_df <- rbind(temp_results_df, adder)
  }
  temp_results_df<-temp_results_df[-1,]
  temp_results_df <- temp_results_df[order(temp_results_df$MSE, decreasing=F),]
  bestMod <- temp_results_df$num[1]
  
  results_df <- rbind(results_df, temp_results_df)
  
  modelList[[i]] <- trainModelList[[bestMod]]
  modelList[[i]]$alpha <- alpha_vals[bestMod]
  
  yhat <- matrix(predict(trainModelList[[bestMod]], s=trainModelList[[bestMod]]$lambda.1se, newx=test_x), 
                 nrow(test_x),ncol(test_y))
  
  for (j in 1:ncol(test_y)){
    condition<-colnames(test_y)[j]
    cc=cor(yhat[,j], test_y[,j])
    iteration=i
    adder <- data.frame(condition = condition, cor = cc, iteration = iteration)
    cor_res_test <- rbind(cor_res_test, adder)
  }
  
  yhat <- matrix(predict(trainModelList[[bestMod]], s=trainModelList[[bestMod]]$lambda.1se, newx=train_x), 
                 nrow(train_x),ncol(train_y))
  for (j in 1:ncol(train_y)){
    condition<-colnames(train_y)[j]
    cc=cor(yhat[,j], train_y[,j])
    iteration=i
    adder <- data.frame(condition = condition, cor = cc, iteration = iteration)
    cor_res_train <- rbind(cor_res_train, adder)
  }
}
cor_res_test <- cor_res_test[-1,]
cor_res_train <- cor_res_train[-1,]

g1 <- ggplot(cor_res_test, aes(reorder(condition, cor, FUN = median), cor, color=condition)) +
  geom_boxplot() + geom_jitter(width=0.1) + ggtitle("Test Correlations") +
  theme(axis.text.x = element_text(size=10, angle = 90, hjust=1,vjust=1), legend.position = 'none')

g2 <- ggplot(cor_res_train, aes(reorder(condition, cor, FUN = median), cor, color=condition)) +
  geom_boxplot() + geom_jitter(width=0.1) + ggtitle("Train Correlations") +
  theme(axis.text.x = element_text(size=10, angle = 90, hjust=1,vjust=1), legend.position = 'none')


pdf("final_corr.pdf", height = 15, width = 10)
grid.arrows(g2,g1,nrow=2)
dev.off()

results_df <- results_df[-1,]
gloss <- ggplot(results_df, aes(x=alpha, y=lambda, color = MSE))+
  geom_point(size=3)

pdf("total_loss.pdf", height = 20, width = 10)
gloss+facet_grid(~iteration, rows = 5)
dev.off()

write.table(results_df, "total_loss.txt", sep='\t', quote=F, row.names=F)

write.table(cor_res_train, "train_corr_table.txt", sep='\t', quote=F, row.names=F)
write.table(cor_res_test, "test_corr_table.txt", sep='\t', quote=F, row.names=F)
save(modelList, file="modelList.RData")