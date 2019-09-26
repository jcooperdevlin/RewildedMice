### Integration Model Test1
### 3.27

setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/IntegrationModel")

library(ggplot2)
#library(MultivariateRandomForest)
library(matrixStats)
library(fastR)
library(ComplexHeatmap)
library(circlize)
library(psych)

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



### are predictors correlated?
overlaps <- intersect(rownames(trainY), rownames(trainXs))
overlaps <- sample(overlaps, length(overlaps))

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
cor_res <- data.frame(condition=NA, cor=NA, iteration=NA)
modelList <- list()
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
  results_df <- data.frame(alpha=NA, lambda=NA,MSE=NA, num=NA)
  alpha_vals <- seq(0,1, length=50)
  trainModelList=list()
  for (k in 1:50) {
    trainModelList[[k]]<- cv.glmnet(train_x, train_y, type.measure="mse",nfolds = 7,
                                    alpha=alpha_vals[k],family="mgaussian", 
                                    standardize.response=T)
    yhat <- matrix(predict(trainModelList[[k]], s=trainModelList[[k]]$lambda.1se, newx=test_x), 
                   nrow(test_x),ncol(test_y))
    mse <- mean((test_y - yhat)^2)
    adder <- data.frame(alpha = round(alpha_vals[k],3), lambda = trainModelList[[k]]$lambda.1se, 
                        MSE = mse, num=k)
    results_df <- rbind(results_df, adder)
  }
  results_df<-results_df[-1,]
  results_df <- results_df[order(results_df$MSE, decreasing=F),]
  bestMod <- results_df$num[1]
  
  modelList[[i]] <- trainModelList[[bestMod]]
  
  yhat <- matrix(predict(trainModelList[[bestMod]], s=trainModelList[[bestMod]]$lambda.1se, newx=test_x), 
                 nrow(test_x),ncol(test_y))
  
  for (j in 1:ncol(test_y)){
    condition<-colnames(test_y)[j]
    cc=cor(yhat[,j], test_y[,j])
    iteration=i
    adder <- data.frame(condition = condition, cor = cc, iteration = iteration)
    cor_res <- rbind(cor_res, adder)
  }
}
cor_res <- cor_res[-1,]

g <- ggplot(cor_res, aes(condition, cor, color=condition)) +
  geom_boxplot() + geom_jitter(width=0.1) +
  theme(axis.text.x = element_text(size=10, angle = 90))

pdf("final_corr.pdf", height = 10, width = 10)
g
dev.off()

write.table(cor_res, "final_corr_table.txt", sep='\t', quote=F, row.names=F)
save(modelList, file="modelList.RData")


#
#

#### ^^^^ ran on cluster...

cor_res <- read.table("final_corr_table.txt", T, '\t')
g <- ggplot(cor_res, aes(reorder(condition, cor, FUN = median), cor, color=condition)) +
  geom_boxplot() + geom_jitter(width=0.1) +
  theme(axis.text.x = element_text(size=10, angle = 90), legend.position = 'none')

load("modelList.RData")

#
#
### check models


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
### garbage
# Prediction=build_forest_predict(trainXs[1:65,], trainY[1:65,], 50, 10, 5, trainXs[71:84,])
# 
# cor(c(Prediction[,10]), c(trainY[61:84,10]))
# 
# #
# #
# # export for sklearn
# write.table(trainXs[1:65,], "train_x.txt", row.names=F, quote=F, sep='\t')
# write.table(trainY[1:65,], "train_y.txt", row.names=F, quote=F, sep='\t')
# write.table(trainXs[66:84,], "test_x.txt", row.names=F, quote=F, sep='\t')
# write.table(trainY[66:84,], "test_y.txt", row.names=F, quote=F, sep='\t')
# 
# 
# 
# 
# 
# #
# #
# #
# # try it by hand to understand
# trainX = trainXs[1:60,]
# trainY = trainY[1:60,]
# n_tree = 10
# m_feature = 10
# min_leaf = 5
# testX = trainXs[61:84,]
# 
# build_forest_predict <- function(trainX, trainY, n_tree, m_feature, min_leaf, testX){
#   if (class(n_tree)=="character" || n_tree%%1!=0 || n_tree<1) stop('Number of trees in the forest can not be fractional or negative integer or string')
#   if (class(m_feature)=="character" || m_feature%%1!=0 || m_feature<1) stop('Number of randomly selected features considered for a split can not be fractional or negative integer or string')
#   if (class(min_leaf)=="character" || min_leaf%%1!=0 || min_leaf<1 || min_leaf>nrow(trainX)) stop('Minimum leaf number can not be fractional or negative integer or string or greater than number of samples')
#   
#   theta <- function(trainX){trainX}
#   results <- bootstrap::bootstrap(1:nrow(trainX),n_tree,theta) 
#   b=results$thetastar
#   
#   Variable_number=ncol(trainY)
#   if (Variable_number>1){
#     Command=2
#   }else if(Variable_number==1){
#     Command=1
#   } 
#   
#   Y_HAT=matrix(  0*(1:Variable_number*nrow(testX)),  ncol=Variable_number,   nrow=nrow(testX)  )
#   Y_pred=NULL
#   
#   for (i in 1:n_tree){
#     Single_Model=NULL
#     X=trainX[ b[ ,i],  ]
#     Y=matrix(trainY[ b[ ,i],  ],ncol=Variable_number)
#     Inv_Cov_Y = solve(cov(Y)) # calculate the V inverse
#     if (Command==1){
#       Inv_Cov_Y=matrix(rep(0,4),ncol=2)
#     }
#     Single_Model=build_single_tree(X, Y, m_feature, min_leaf,Inv_Cov_Y,Command)
#     Y_pred=single_tree_prediction(Single_Model,testX,Variable_number)
#     for (j in 1:Variable_number){
#       Y_HAT[,j]=Y_HAT[,j]+Y_pred[,j]
#     }
#   }
#   Y_HAT=Y_HAT/n_tree
#   return(Y_HAT)
# }
# 
# 
# 
# 
# 
# 
# #
# #
# #
# ## manual
# 
# #Input and Output Feature Matrix of random data (created using runif)
# trainX=matrix(runif(70*100),70,100)
# trainY=matrix(runif(70*5),70,5)
# n_tree=500
# m_feature=2
# min_leaf=1
# #testX=matrix(runif(10*100),10,100)
# #Prediction size is 10 x 5, where 10 is the number
# #of testing samples and 5 is the number of output features
# Prediction=build_forest_predict(trainX[1:50,], trainY[1:50,], n_tree, m_feature, min_leaf, trainX[51:70,])
# Compare=trainY[51:70,]
# 
# cor(c(Prediction[,4]), c(Compare[,4]))
# 
# 
#
#
#
#
#
#
### glmnet usage

# Gaussian
x=matrix(rnorm(100*20),100,20)
y=rnorm(100)
fit1=glmnet(x,y)
print(fit1)
coef(fit1,s=0.01) # extract coefficients at a single value of lambda
predict(fit1,newx=x[1:10,],s=c(0.01,0.005)) # make predictions
#multivariate gaussian
y=matrix(rnorm(100*3),100,3)
fit1m=glmnet(x,y,family="mgaussian")
plot(fit1m,type.coef="2norm")

fit1=glmnet(x=trainXs[1:65,], y=trainY[1:65,], family="mgaussian")
aa=fit1$a0
preds=predict(fit1, trainXs[66:84,])
preds1 <- matrix(preds, 19*100, 10)
