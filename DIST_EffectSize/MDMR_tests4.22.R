#### MDMR Tests for Rewilding Data
#### Multivariate Distance Matrix Regression

library(MDMR)

### example

# Source data
data(mdmrdata)
# Get distance matrix
D <- dist(Y.mdmr, method = 'euclidean')
# Conduct MDMR
mdmr.res <- mdmr(X = X.mdmr, D = D)
summary(mdmr.res)
################################################################
## Conducting MDMR on data comprised of dependent observations
################################################################
# Source data
data("clustmdmrdata")
# Get distance matrix
D <- dist(Y.clust)
# Conduct mixed-MDMR
mixed.res <- mixed.mdmr(~ x1 + x2 + (x1 + x2 | grp),
                        data = X.clust, D = D)
summary(mixed.res)

#
#
#
#
#
data(mdmrdata)
# --- Method 1 --- #
es_res = delta(X.mdmr, Y = Y.mdmr, dtype = "euclidean", niter = 1, seed = 12345, plot.res = T)
# --- Method 2 --- #
D <- dist(Y.mdmr, method = "euclidean")
G <- gower(D)
q <- ncol(Y.mdmr)
G.list <- vector(mode = "list", length = q)
names(G.list) <- names(Y.mdmr)
for(i in 1:q) {
  Y.shuf <- Y.mdmr
  Y.shuf[,i] <- sample(Y.shuf[,i])
  G.list[[i]] <- gower(dist(Y.shuf, method = "euclidean"))
}
delta(X.mdmr, G = G, G.list = G.list, plot.res = T)

res <- pcoa(D)
plotter <- data.frame(res$vectors, col=X.mdmr[,3])

ggplot(plotter, aes(Axis.1, Axis.2, color=col)) + 
  geom_point()



#
#

#
#
## let's try with our data

mice_data <- read.table("data/FACS_data/BLOOD_lymph_FACS_metadata_names.txt", sep='\t', header=T)
name_change <- read.table("data/FACS_data/lymph_name_change.txt")
colnames(mice_data)[13:27] <- as.character(name_change$V2)

pc_use <- mice_data[,13:27]

D <- dist(pc_use, method = "euclidean")
res <- pcoa(D)
plotter <- data.frame(res$vectors, Environment=mice_data$Environment, Genotype = mice_data$Genotype)

ggplot(plotter, aes(Axis.1, Axis.2, color=Environment)) + 
  geom_point()

effectors <- mice_data[,c(4,5,6,8,9,10)]
mdmr.res <- mdmr(X = effectors, D = D)
summary(mdmr.res)

delta(effectors, Y = pc_use, dtype = "euclidean", niter = 1, seed = 12345, plot.res = T)



#
#
#
cyt_data <- readRDS("IntegrationModel/ml_inputs/cytokines.RData")
meta <- readRDS("IntegrationModel/ml_inputs/mouse_metadata.RData")
meta$Env_gen <- paste0(meta$Environment, ":", meta$Genotype)

D <- dist(cyt_data, method = "euclidean")
res <- pcoa(D)
plotter <- data.frame(res$vectors, Environment=meta$Environment, 
                      Genotype = meta$Genotype, combo = meta$Env_gen)

ggplot(plotter, aes(Axis.1, Axis.2, color=Environment, shape=Genotype)) + 
  geom_point()

effectors <- meta[,c(2,3,4,6,7,8)]
mdmr.res <- mdmr(X = effectors, D = D)
summary(mdmr.res)

delta_res=delta(effectors, Y = cyt_data, dtype = "euclidean", niter = 1, seed = 12345, plot.res = T)
delta_res <- data.frame(t(delta_res))
