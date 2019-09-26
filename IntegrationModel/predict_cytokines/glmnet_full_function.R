##### a function to run GLMNet with certain X and certain Y
library(caret)
library(glmnet)
library(ggplot2)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(psych)
library(gridExtra)
#library(reshape)

varImp <- function(object, lambda = NULL, ...) {
  
  ## skipping a few lines
  
  beta <- predict(object, s = lambda, type = "coef")
  if(is.list(beta)) {
    out <- do.call("cbind", lapply(beta, function(x) x[,1]))
    out <- as.data.frame(out)
  } else out <- data.frame(Overall = beta[,1])
  out <- abs(out[rownames(out) != "(Intercept)",,drop = FALSE])
  out
}

run_glm <- function(x, y, output_dir, type, p_threshold=0.05){
  if(dir.exists(output_dir)==F){
    dir.create(file.path(output_dir), showWarnings = FALSE)
  }
  
  trainXs <- x
  trainY <- y
  
  overlaps <- intersect(rownames(trainY), rownames(trainXs))
  overlaps <- sample(overlaps, length(overlaps))
  
  trainXs <- trainXs[overlaps,]
  trainY <- trainY[overlaps,]
  
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
  
  cor_mat <- matrix(cor_adder, ncol(trainXs), ncol(trainY))
  rownames(cor_mat) <- colnames(trainXs)
  colnames(cor_mat) <- colnames(trainY)
  
  p_adder <- p.adjust(p_adder, method="BH")
  p_mat <- matrix(p_adder, ncol(trainXs), ncol(trainY))
  rownames(p_mat) <- colnames(trainXs)
  colnames(p_mat) <- colnames(trainY)
  
  
  good_genes <- rownames(cor_mat)[unique(which(p_mat <= p_threshold, arr.ind = T)[,1])]
  
  if(length(good_genes)>3){
    namer <- paste0(output_dir, "/", type, "_feature_vector.txt")
    good_gene_tab <- data.frame(num=1:length(good_genes), Feature=good_genes)
    write.table(good_gene_tab, namer, row.names = F, sep='\t', quote=F)
    
    namer <- paste0(output_dir, "/", type, "_cyt_corr_heat.png")
    png(namer, height =5, width = 10, units = 'in', res=150)
    hh=Heatmap(cor_mat[good_genes,], show_row_names = F, show_column_names = F,
               col = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")))
    print(hh)
    dev.off()
    
    print(paste0(length(good_genes), " ", type, " selected to train GLM"))
    #
    #
    
    ### RANDOMLY TEST 10 Times different test and train sets for full dataset with all genes
    
    results_df <- data.frame(alpha=NA, lambda=NA,MSE=NA, num=NA)
    alpha_vals <- seq(0.01,1, length=20)
    #lambda_vals <- seq(0.001,0.1,by = 0.001)
    lambda_vals <- seq(0.0001, 5, length = 100)
    trainModelList=list()
    for(m in 1:10) {
      if((m) %% 1==0){print(paste0("Starting iteration...", m, " of 10"))}
      for (k in 1:20) {
        if((k) %% 5==0){print(paste0("Starting subiteration...", k, " of 20"))}
        modeler<- cv.glmnet(trainXs, trainY, type.measure="mse",nfolds = 7,
                            alpha=alpha_vals[k],family="mgaussian", lambda = lambda_vals,
                            standardize.response=T)
        
        adder <- data.frame(alpha = rep(round(alpha_vals[k],3),length(modeler$lambda)), lambda = modeler$lambda, 
                            MSE = modeler$cvm, num=rep(m,length(modeler$lambda)))
        results_df <- rbind(results_df, adder)
      }
    }
    results_df<-results_df[-1,]
    results_df <- results_df[order(results_df$MSE, decreasing=F),]
    
    results_df$Log2_MSE <- log2(results_df$MSE+1)
    results_df$combo <- as.factor(paste0(results_df$alpha, ";", round(results_df$lambda,3)))
    gloss <- ggplot(results_df, aes(x=as.factor(alpha), y=Log2_MSE, color = lambda))+
      geom_boxplot(outlier.shape = NA) + geom_jitter() + 
      stat_summary(fun.y = 'median', geom = 'point', color = "dodgerblue2", shape=3) +
      ylab("Log2 MSE") + xlab("Alpha") +
      scale_color_gradient2(low="grey90", mid="dodgerblue", high = "navyblue", 
                            midpoint=median(results_df$lambda)) +
      #annotate('point', x=best_alpha, y=log2(best_lambda+1), color="red3", size=3) +
      theme_bw()
    
    agg_res <- aggregate(.~combo, data=results_df, mean)
    agg_res <- agg_res[order(agg_res$MSE, decreasing=F),]
    best_alpha <- agg_res$alpha[1]
    best_lambda <- agg_res$lambda[1]
    #plot_combo <- agg_res$combo[100]
    lowest_MSE <- agg_res$Log2_MSE[1]
    
    
    namer <- paste0(output_dir, "/", type, "_total_loss.pdf")
    pdf(namer, height = 5, width = 10)
    print(
      gloss + annotate('text', x="0.2", 
                       y=quantile(results_df$Log2_MSE, 0.99),
                       label = paste0("best alpha = ", best_alpha, '\n',
                                      "best lambda = ", round(best_lambda,3)),
                       color="red3", size=5) +
        annotate('point', x=as.character(best_alpha), 
                 y=lowest_MSE, color = "red3", size=3)
    )
    dev.off()
    
    namer <- paste0(output_dir, "/", type, "_total_loss.txt")
    write.table(results_df, namer, row.names=F, sep='\t', quote=F)
    
    
    print(paste0("Alpha selected at...", best_alpha))
    print(paste0("Lambda selected at...", round(best_lambda,3)))
    print(paste0(""))
    
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
      test_x <- trainXs[71:length(overlaps),]
      test_y <- trainY[71:length(overlaps),]
      
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
    
    namer <- paste0(output_dir, "/", type, "_final_corr.pdf")
    pdf(namer, height = 6, width = 14)
    print(g1)
    dev.off()
    
    namer <- paste0(output_dir, "/", type, "_final_corr.txt")
    write.table(cor_res, namer, sep='\t', quote=F, row.names=F)
    namer <- paste0(output_dir, "/", type, "_modelList.Rdata")
    save(modelList, file=namer)
    
    avg_importance <- data.frame(varImp(modelList[[1]], lambda=modelList[[1]]$lambda)[1,],
                                 Feature=NA)
    for(n in 1:length(modelList)){
      varImp_df <- varImp(modelList[[n]], lambda=modelList[[n]]$lambda)
      varImp_df$Feature <- rownames(varImp_df)
      avg_importance <- rbind(avg_importance, varImp_df)
    }
    avg_importance <- avg_importance[-1,]
    avg_mat <- aggregate(.~Feature, data=avg_importance, mean)
    
    if(max(avg_mat[,-1])>0){
      
      avg_mat2 <- avg_mat[which(rowSums(avg_mat[,-1])>0.01),]
      
      if(nrow(avg_mat2)>0){
        avg_mat <- avg_mat2
      } else {
        avg_mat <- avg_mat
      }
      
      avg_heat <- avg_mat[,-1]
      rownames(avg_heat) <- avg_mat[,1]
      namer <- paste0(output_dir, "/", type, "_feature_matrix.txt")
      write.table(avg_mat, namer, row.names=F, sep='\t', quote=F)
      
      #fix funky rownames
      for(b in 1:nrow(avg_heat)){
        if(grepl(';',rownames(avg_heat)[b])){
          lookerupper <- paste0("^",rownames(avg_heat)[b], "$")
          if(grepl("\\[", lookerupper)){
            lookerupper1 <- strsplit(lookerupper, "\\[")
            lookerupper2 <- strsplit(lookerupper, "\\]")
            lookerupper <- paste0(lookerupper1[[1]][1], "\\[.*\\]",
                                  lookerupper2[[1]][2])
          }
          lookupnum <- grep(lookerupper, good_gene_tab$Feature, fixed=F)
          rownames(avg_heat)[b] <- paste0("OTU_", lookupnum)
        }
      }
      
      
      labs = colnames(avg_heat)
      ll <- columnAnnotation(labels = anno_text(labs, gp=gpar(fontsize=8),
                                                which = "column", rot=60, 
                                                just = 'right', offset=unit(3.5,"cm")), 
                             height = unit(2,"cm"))
      
      heighter <- c(0.25*nrow(avg_heat))+3.5+2
      namer <- paste0(output_dir, "/", type, "_feat_imp_heat.png")
      png(namer, height = heighter, width = 20, units = 'cm', res=150)
      hh=Heatmap(avg_heat, show_row_names = T, show_column_names = F,
                 name="Importance",
                 col = colorRamp2(c(0, max(avg_heat)), c("white", "navyblue")),
                 row_names_gp = gpar(fontsize=8),
                 bottom_annotation = ll,
                 bottom_annotation_height = unit(3.5,"cm"))
      print(hh)
      dev.off()
    } else {
      print("All Feature importances were 0...exiting")
    }
  } else {
    print("Less than 3 features were significantly correlated...exiting")
  }
}