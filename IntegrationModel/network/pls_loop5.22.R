### pls loop

library(mixOmics)

genes = readRDS("ml_inputs/gene_modules.RData")
cytokines = readRDS("ml_inputs/cytokines.RData")
cells = readRDS("ml_inputs/facs.RData")
otus = readRDS("ml_inputs/otus_relab.RData")

inputs = list(genes=genes, cytokines=cytokines, cells=cells, otus=otus)
combos = combn(4,2)

for ( i in 1:ncol(combos)){
  X=inputs[[combos[1,i]]]
  X_name=names(inputs)[combos[1,i]]
  Y=inputs[[combos[2,i]]]
  Y_name=names(inputs)[combos[2,i]]
  
  first.spls <- spls(X, Y, ncomp = 10)
  first.val <- perf(first.spls, validation = "Mfold", folds = 5, repeats = 10)
  ncomp = which(first.val$Q2 < 0.0975)[1]
  
  true.spls <- spls(X, Y, ncomp = ncomp)
  
  cord.X = cor(true.spls$X, true.spls$variates$X[, 1:ncomp], use = "pairwise")
  cord.Y = cor(true.spls$Y, true.spls$variates$Y[, 1:ncomp], use = "pairwise")
  true_sim = cord.X %*% t(cord.Y)
  
  true_sim <- data.frame(item1=rownames(true_sim), true_sim)
  namer <- paste0("sim_mats/", X_name, "_vs_", Y_name, ".txt")
  write.table(true_sim, namer, sep='\t', row.names=F, quote=F)
  
  true_flat <- c(cord.X %*% t(cord.Y))
  true_flat[is.na(true_flat)] <- 0
  
  num_sim <- 10
  null_mat <- matrix(NA, 1, ncol(X)*ncol(Y))
  for(j in 1:num_sim){
    null.spls <- spls(X, Y[sample(rownames(Y)),], ncomp = ncomp)
    cord.X = cor(null.spls$X, null.spls$variates$X[, 1:ncomp], use = "pairwise")
    cord.Y = cor(null.spls$Y, null.spls$variates$Y[, 1:ncomp], use = "pairwise")
    null_sim = cord.X %*% t(cord.Y)
    
    null_flat <- c(null_sim)
    null_mat <- rbind(null_mat, null_flat)
  }
  null_mat <- null_mat[-1,]
  null_mat[is.na(null_mat)] <- 0
  
  p.value = NULL
  for (k in 1:ncol(null_mat)){
    p=mean(abs(null_mat[,k])>=abs(true_flat[k]))
    p.value = c(p.value,p)
  }
  p.adj = p.adjust(p.value, "BH")
  p_mat <- matrix(p.adj, ncol(X), ncol(Y))
  rownames(p_mat) <- colnames(X)
  colnames(p_mat) <- colnames(Y)
  
  p_mat <- data.frame(item1=rownames(p_mat), p_mat)
  namer <- paste0("sim_mats/", X_name, "_vs_", Y_name, "_p.adj.txt")
  write.table(p_mat, namer, sep='\t', row.names=F, quote=F)
}