#### Try to build a PLS network of our data types

#### gene modules, OTUs, Immune cells, Cytokines

setwd("/Volumes/lokep01lab/lokep01labspace/Rewilding_Data/int/IntegrationModel/network/")

library(ggplot2)
library(igraph)
#library(mixOmics)
library(jmuOutlier)



spls_ex <- readRDS("spls_ex.RData")
huh=spls_ex$variates[[1]]
huh=spls_ex$explained_variance

huh=spls_ex$variates[[1]]
huh2=spls_ex$variates[[2]]



net_net = readRDS("net-1.RData")
corr_mat = net_net$M_mrna_mirna

mat = spls_ex
comp=c(1:2)

cord.X = cor(mat$X$mrna, mat$variates$mrna[, comp], use = "pairwise")
cord.Y = cor(mat$X$mirna, mat$variates$mirna[, comp], use = "pairwise")
corr_mat = cord.X %*% t(cord.Y)
corr_flat <- c(corr_mat)

num_sim=50
sim_test = matrix(0, 1, 36800)
for(i in 1:num_sim){
og_x <- mat$X$mrna
sim_x <- og_x[sample(rownames(og_x)),]
cord.X = cor(mat$X$mrna, mat$variates$mrna[, comp], use = "pairwise")
cord.X = corr.test(mat$X$mrna, mat$variates$mrna[, comp])

rs = cord.X$r
ps = cord.X$p

og_y <- mat$X$mirna
sim_y <- og_y[sample(rownames(og_y)),]
cord.Y = cor(sim_y, mat$variates$mirna[, comp], use = "pairwise")
cord.Y = corr.test(mat$X$mirna, mat$variates$mirna[, comp])

rs2 = cord.Y$r
ps2 = cord.Y$p

sim_mat = cord.X %*% t(cord.Y)
sim_test = rbind(sim_test, c(sim_mat))
}
p_test=ps %*% t(ps2) 

p.value = NULL
for (j in 1:ncol(sim_test)){
  p=mean(abs(sim_test[,j])>=abs(corr_flat[j]))
  p.value = c(p.value,p)
}
summary(p.adjust(p.value, "BH"))

plot(density(p.adjust(p.value, "BH")))
abline(v=0.05, col='red')

#
#
#
#

#
#
#
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

#PLS
result <- pls(X, Y, ncomp = 3)  # where ncomp is the number of dimensions/components to choose
tune.pls <- perf(result, validation = 'loo', criterion = 'all', progressBar = FALSE)

# SPLS
ncomp = 10
result.spls <- spls(X, Y, ncomp = ncomp, keepX = c(rep(10, ncomp)), mode = 'regression')
tune.spls <- perf(result.spls, validation = 'Mfold', folds = 10,
                  criterion = 'all', progressBar = FALSE)